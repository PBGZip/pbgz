
#include "compress_engine.h"
#include "fastq_actuator.h"
#include "binary_actuator.h"
#include "utils/path_util.h"
#include "coder_ppmd.h"
#include "coder_json.h"

int32_t CompressEngine::init() {
    if (0 != PbgzEngine::init()) {
        return -1;
    }

    indexBlockQueue.setCapility(2);
    freeIndexBlockQueue.setCapility(2);
    for (uint32_t i = 0; i < freeIndexBlockQueue.getCapility(); ++i) {
        RoughIOBlock* ptr = MemoryUtil::safeNewClass<RoughIOBlock>();
        if (ptr == nullptr) {
            LOG_ERROR("Failed to create RoughIOBlock for freeIndexBlockQueue");
            return -1;
        }
        freeIndexBlockQueue.push(ptr);
    }
    
    return 0;
}

CompressEngine::~CompressEngine() {
    while(!freeIndexBlockQueue.empty()) {
        RoughIOBlock* ptr = freeIndexBlockQueue.get();
        MemoryUtil::safeDeleteClass(ptr);
    }

    while(!indexBlockQueue.empty()) {
        RoughIOBlock* ptr = indexBlockQueue.get();
        MemoryUtil::safeDeleteClass(ptr);
    }
}

BlockReader* CompressEngine::createBlockReader() {
    BlockReader* blockReader = MemoryUtil::safeNewClass<BlockReader>(ioReader);
    if (blockReader == nullptr) {
        LOG_ERROR("Create block reader failed.");
        return nullptr;
    }
    if (0 != blockReader->init()) {
        LOG_ERROR("BlockReader init failed");
        MemoryUtil::safeDeleteClass(blockReader);
        return nullptr;
    }

    return blockReader;
}

BlockWriter* CompressEngine::createBlockWriter() {
   PbgzBlockWriter* blockWriter = MemoryUtil::safeNewClass<PbgzBlockWriter>(ioWriter);
    if (blockWriter == nullptr) {
        LOG_ERROR("pbgzWriter is NULL.");
        return nullptr;
    }
    blockWriter->init();
    blockWriter->setBaseFileMeta(baseFileMeta);
    blockWriter->writeBaseFileMeta();
    return blockWriter;
}

Actuator* CompressEngine::actuatorPreProc(Actuator* actuator, RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
    FastqActuator* fastqActuator = dynamic_cast<FastqActuator*>(actuator);
    if (fastqActuator != nullptr) {
        if (0 != fastqActuator->preAnalysis()) {
            LOG_INFO("Fastq preAnalysis failed");
            MemoryUtil::safeDeleteClass(fastqActuator);
            return MemoryUtil::safeNewClass<BinaryActuator>(inBlockPtr, outBlockPtr);
        }
    }
    return actuator;
}

int32_t CompressEngine::actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock* outBlockPtr) {
    int32_t ret = actuator->compress();
    return ret;
}


int32_t CompressEngine::engineStartPreProc() {
    if (!initReference()) {
        LOG_ERROR("Engine init reference failed.");
        return -1;
    }
    return 0;
}

bool CompressEngine::initReference() {
    if (!parameter.referenceGenic.empty()){
        pRefGene = MemoryUtil::safeNewClass<Reference>(parameter.referenceGenic, parameter.threadNum);
    }
    if (pRefGene) {
        if (!pRefGene->makeIndex()) {
            MemoryUtil::safeDeleteClass(pRefGene);
        }
    }
    // Refresh file meta
    if (pRefGene) {
        Json::Value refeMeta;
        refeMeta["squash_len"] = (Json::Value::Int64)(pRefGene->getSquashLength());
        refeMeta["fasta_name"] = PathUtil::getFileName(pRefGene->getFastaFileName());
        refeMeta["fasta_len"] = (Json::Value::Int64)(PathUtil::getFileSize(pRefGene->getFastaFileName()));
        refeMeta["fasta_md5"] = pRefGene->getFastaChecksum();
        refeMeta["ni_name"] = PathUtil::getFileName(pRefGene->getNiFilePath()); /* Contains md5 information, used for decompression verification */
        refeMeta["max_block_len"] = 16 << 20;
        bool isPackRefeInHeader = !parameter.isUnpackRef && parameter.outputFile == "/dev/stdout";
        if (isPackRefeInHeader) {
            refeMeta["blocks"] = calcPackRefeBlockSize();
            baseFileMeta.setMetaData("refe",refeMeta);
            int64_t maxRefLen = 0;
            int64_t totolEncLen = 0;
            (void)packReference(maxRefLen, totolEncLen, false);
        } else {
            refeMeta["blocks"] = 0;
            baseFileMeta.setMetaData("refe",refeMeta);
        }
    }
    return true;
}

uint32_t CompressEngine::calcPackRefeBlockSize() {
    int64_t each = (16 << 20);  
    int64_t total = pRefGene->getSquashLength();

    uint32_t result = total / each;
    int64_t remain = total % each;
    if (remain > 0) {
        result++;
    }
    return result;
}


/*  Save reference genome */
int64_t CompressEngine::packReference(int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef) {
    int64_t block = 0, offset = 0;
    std::vector<std::thread> tpools;
    
    uint32_t pcnt = parameter.threadNum;
    int64_t each = (16 << 20);
    int64_t total = pRefGene->getSquashLength();
    int64_t remain = total;
    maxBlockLen =0;
    totalEncLen = 0;
    
    BlockingQueue<RefeInfo> inputPool(pcnt);
    Reference* refe = pRefGene;
    int64_t current;
    std::mutex m;
    BlockingQueue<RoughIOBlock*>& freeOutput = freeOutputPool;
    BlockingQueue<RoughIOBlock*>& outputPool = outputDataPool;
    uint32_t count = blockCount;

    // Start threads
    for (uint32_t idx = 0; idx < pcnt; ++idx) {
        tpools.push_back(std::thread([&inputPool, &m, &freeOutput, &outputPool, &refe, &each, &maxBlockLen, &totalEncLen, &isSanitizeRef, &count]() {
            uint8_t* output = MemoryUtil::safeAlloc<uint8_t>(each);
            while(true) {
                RefeInfo refe2do = inputPool.get();
                int64_t plen = refe2do.second.second;
                if (plen == 0) {
                    break;
                }
                const uint8_t *p = refe->getSquash() + refe2do.second.first;
                if (isSanitizeRef) {
                    refe->sanitizeRefSquash(refe2do.second.first, plen);
                }

                coder_io refeIo(output, each);
                coder_ppmd cppmd(&refeIo);
                cppmd.encode(p, plen);
                cppmd.encode_flush();

                Json::Value meta;
                meta["block"] = (Json::Value::Int)refe2do.first.first;
                meta["offset"] = (Json::Value::Int64)refe2do.second.first;
                meta["srclen"] = (Json::Value::Int)plen;
                meta["dstlen"] = (Json::Value::Int)refeIo.data_len;
                meta["coder"] = refeIo.meta;
                coder_json cmeta;
                std::string metaString;
                cmeta.encoder(meta, metaString); /*  Compress block meta */
            
                int64_t currBlockLen = metaString.length() + refeIo.data_len;
                if (currBlockLen >= plen) {
                    LOG_ERROR("reference block pack failed");
                    break;
                }
                /* Write compressed stream of current block's meta */
                RoughIOBlock* outBlock = freeOutput.get();
                outBlock->reset();
                memcpy(outBlock->getCurrent(), refeIo.data, refeIo.data_len);
                outBlock->setDataLen(refeIo.data_len);
                memcpy(outBlock->getCurrent(), metaString.c_str(), metaString.length());
                outBlock->setMetaLen(metaString.length());
                outBlock->setBlockType(REFERENCE);
                outBlock->setBlockId(count + refe2do.first.first);
                outputPool.push(outBlock);

                m.lock();
                maxBlockLen = (currBlockLen > maxBlockLen) ? currBlockLen : maxBlockLen;
                totalEncLen += currBlockLen;
                m.unlock();
            }
            MemoryUtil::safeFree(output);
        }));
    }

    while (remain > 0) {
        current = std::min(each, remain);
        inputPool.push(std::make_pair(std::make_pair(block, nullptr), std::make_pair(offset, current)));
        block++;
        offset += current;
        remain -= current;
    }
    // Write end marker
    for (uint32_t i = 0; i < pcnt; ++i) {
        std::pair<std::pair<int64_t, uint8_t *>, std::pair<int64_t, int64_t>> refe2do;
        refe2do.second.second = 0;
        inputPool.push(refe2do);
    }
    for (auto &t : tpools) {
        if (t.joinable()) {
            t.join();
        }
    }

    return block;
}

int32_t CompressEngine::engineStartAfterProc() {
    bool isPackRefeInTail = pRefGene && !parameter.isUnpackRef && parameter.outputFile != "/dev/stdout";
    if (isPackRefeInTail) {
        int64_t maxRefLen = 0;
        int64_t totolEncLen = 0;
        Json::Value refeMeta;
        refeMeta["squash_len"] = (Json::Value::Int64)(pRefGene->getSquashLength());
        refeMeta["fasta_name"] = PathUtil::getFileName(pRefGene->getFastaFileName());
        refeMeta["fasta_len"] = (Json::Value::Int64)(PathUtil::getFileSize(pRefGene->getFastaFileName()));
        refeMeta["fasta_md5"] = pRefGene->getFastaChecksum();
        refeMeta["ni_name"] = PathUtil::getFileName(pRefGene->getNiFilePath());
        refeMeta["offset"]= 0;
        dynamicFileMeta.setMetaData("refe", refeMeta);
        FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
        if (fileWriter != nullptr) {
            int64_t refBlockCount = packReference(maxRefLen, totolEncLen);
            dynamicFileMeta.getMetaData("refe")["max_block_len"] = maxRefLen;
            dynamicFileMeta.getMetaData("refe")["blocks"] = refBlockCount;
        }
    }
    return 0;
}

void CompressEngine::setDataBlockPosition(uint32_t blockId) {
    if (!parameter.isMakeIndex) {
        return;
    }
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    if (fileWriter == nullptr) {
        return;
    }
    pbgzIndex.setBlockPosition(blockId, fileWriter->getCurrentPos());
    return;
 }

 void CompressEngine::startWriteIndexTask() {
    if (!parameter.isMakeIndex) {
        return;
    }

    return;
 }
