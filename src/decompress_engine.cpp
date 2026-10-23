/*
 * decompress_engine.cpp - Source file for pbgz project
 * Copyright (C) 2025 PBGZip
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <set>

#include "decompress_engine.h"
#include "pbgz_index.h"   /* SamIndex, used by the index paths below */
#include "utils/path_util.h"
#include "coder_json.h"
#include "coder_ppmd.h"
#include "pbgz_types.h"
#include "codec_actuator_adapter.h"
#include "sam_info.h"

#include <bzlib.h>

BlockReader* DecompressEngine::createBlockReader() {
    PbgzBlockReader* pbgzReader = MemoryUtil::safeNewClass<PbgzBlockReader>(ioReader);
    pbgzReader->init();
    baseFileMeta = pbgzReader->getBaseFileMeta();
    // Block-size upper bound written during compression, for the actuator's output buffer
    if (baseFileMeta.getMetaData().isMember("block_size")) {
        fileBlockSize = baseFileMeta.getMetaData("block_size").asUInt();
    }
    dynamicFileMeta = pbgzReader->getDynamicFileMeta();
    if (dynamicFileMeta.getMetaData().isMember("qual_prior")) {
        (void)unpackQualPrior(pbgzReader, dynamicFileMeta.getMetaData("qual_prior"));
    }
    /* POS delta prior: hex-encoded in the file meta, parsed once here. */
    if (dynamicFileMeta.getMetaData().isMember("pos_prior")) {
        const std::string& hex = dynamicFileMeta.getMetaData("pos_prior").asString();
        auto nib = [](char c) -> int {
            if (c >= '0' && c <= '9') return c - '0';
            if (c >= 'a' && c <= 'f') return c - 'a' + 10;
            if (c >= 'A' && c <= 'F') return c - 'A' + 10;
            return -1;
        };
        bool ok = (hex.size() % 2 == 0);
        std::vector<uint8_t> prior;
        prior.reserve(hex.size() / 2);
        for (size_t i = 0; ok && i + 1 < hex.size(); i += 2) {
            int h = nib(hex[i]), l = nib(hex[i + 1]);
            if (h < 0 || l < 0) {
                ok = false;
                break;
            }
            prior.push_back((uint8_t)((h << 4) | l));
        }
        if (ok && !prior.empty()) {
            posPriorBlob = std::make_shared<const std::vector<uint8_t>>(std::move(prior));
            LOG_INFO("Loaded POS prior: %zu bytes", posPriorBlob->size());
        }
    }

    if (!initRefGene(pbgzReader)) {
        LOG_INFO("Init reference for decompress failed");
        MemoryUtil::safeDeleteClass(pbgzReader);
        return nullptr;
    }

    do {
        if (parameter.inputFile != STDIN && !parameter.refeGenePos.empty()) {
            if (!initRefeIndex()) {
                LOG_INFO("Init reference index for decompress failed");
                break;
            }

            // Parse refeGenePos (chr1:100-200) into refPosChrName/refPosBegin/refPosEnd
            size_t colonPos = parameter.refeGenePos.find(':');
            if (colonPos == std::string::npos) {
                refPosChrName = parameter.refeGenePos;
                break;
            }
            refPosChrName = parameter.refeGenePos.substr(0, colonPos);

            size_t dashPos = parameter.refeGenePos.find('-');
            if (dashPos == std::string::npos) {
                break;
            }
            refPosBegin = std::stoi(parameter.refeGenePos.substr(colonPos + 1, dashPos - colonPos - 1));
            refPosEnd = std::stoi(parameter.refeGenePos.substr(dashPos + 1));

            LOG_DEBUG("refPosChrName = %s, refPosBegin = %d, refPosEnd = %d", refPosChrName.c_str(), refPosBegin, refPosEnd);
        }
    } while (0);

    return pbgzReader;
}

BlockWriter* DecompressEngine::createBlockWriter() {
    if (parameter.isDecToBam) {
        /* -b: write BAM. BGZF deflate workers are sized at one per three coding
           threads, rounded up (see bamDeflateThreadNum). */
        BlockWriter* blockWriter = MemoryUtil::safeNewClass<BamWriter>(
            ioWriter, bamDeflateThreadNum(parameter.threadNum));
        if (blockWriter == nullptr) {
            LOG_ERROR("Failed to create bam writer.");
            return nullptr;
        }
        return blockWriter;
    }
    BlockWriter* blockWriter = BlockFactory::createBlockWriter(ioWriter);
    if (blockWriter == nullptr) {
        LOG_ERROR("Failed to create block writer.");
        return nullptr;
    }
    return blockWriter;
}

void DecompressEngine::releaseBlockReader(BlockReader* &blockReader) {
    MemoryUtil::safeDeleteClass(blockReader);
}

void DecompressEngine::releaseBlockWriter(BlockWriter* &blockWriter) {
    /* Under -b, finish the BamWriter: flush the residual BGZF block and append the EOF marker. */
    BamWriter* bamWriter = dynamic_cast<BamWriter*>(blockWriter);
    if (bamWriter != nullptr && 0 != bamWriter->finish()) {
        LOG_ERROR("Bam writer finish failed.");
    }
    MemoryUtil::safeDeleteClass(blockWriter);
}

/* Load the QUAL prior on demand for random access: seek to its package-relative offset,
   read the block, and claim it. The read position is restored afterward. */
bool DecompressEngine::unpackQualPrior(PbgzBlockReader* blockReader, Json::Value& priorMeta) {
    if (blockReader == nullptr || !priorMeta.isMember("offset")) {
        return false;
    }

    FileReader* fileReader = dynamic_cast<FileReader*>(ioReader);
    if (fileReader == nullptr) {
        LOG_DEBUG("Input is not seekable, qual prior stays on the streaming path");
        return false;
    }

    const int64_t offset = priorMeta["offset"].asInt64();
    const int64_t packageStart = (int64_t)blockReader->getCurrentFileStart();
    const int64_t blockAddress = packageStart + offset;
    const int64_t packageIndex = blockReader->getCurrentFileIndex();
    /* The offset is only used to seek; registration and lookup use the package index. */
    if (qualPriorConsumer.forPackage(packageIndex) != nullptr) {
        return true;
    }

    const uint32_t dstLen = (uint32_t)priorMeta["dstlen"].asUInt();
    RoughIOBlock* block = MemoryUtil::safeNewClass<RoughIOBlock>(dstLen + (1u << 16));
    if (block == nullptr) {
        return false;
    }

    const size_t readPos = fileReader->getCurrentPos();
    fileReader->seekIO(blockAddress);
    bool ok = (blockReader->readBlock(block) > 0) &&
              qualPriorConsumer.claim(block, packageIndex);
    fileReader->seekIO((int64_t)readPos);
    MemoryUtil::safeDeleteClass(block);

    if (!ok) {
        LOG_ERROR("Read qual prior at offset %ld failed", (long)blockAddress);
    }
    return ok;
}

bool DecompressEngine::initRefGene(PbgzBlockReader* blockReader) {
    if (blockReader == nullptr) {
        return false;
    }

    /* refe metadata may live in baseFileMeta or dynamicFileMeta. */
    bool hasRefe = baseFileMeta.getMetaData("refe").isObject();
    if (!hasRefe) {
        hasRefe = dynamicFileMeta.getMetaData().isMember("refe");
    }
    if (!hasRefe) {
        LOG_DEBUG("No reference");
        return true;
    }

    Json::Value& metaRefe = dynamicFileMeta.getMetaData().isMember("refe")
        ? dynamicFileMeta.getMetaData("refe")
        : baseFileMeta.getMetaData("refe");

    std::string fastaName = metaRefe["fasta_name"].asString();
    int64_t fastaLength = metaRefe["fasta_len"].asInt64();

    LOG_DEBUG("Reference blocks = %d", metaRefe["blocks"].asInt64());

    if (metaRefe["blocks"].asInt64() > 0) {
        /*  There is pack reference in compression package, directly unpack */
        if(unpackReference(blockReader, metaRefe)) {
            return true;
        }
    }

    std::string fastaNameInput = parameter.referenceGenic;
    if (fastaNameInput.empty()) { /* No reference file specified */
        printFastqFileNotMatchInfo(metaRefe);
        LOG_ERROR("reference file needs to be specified to complete decompression");
        return false;
    }

    fastaNameInput = PathUtil::getAbsPath(fastaNameInput);
    int64_t fastqFileLenInput = PathUtil::getFileSize(fastaNameInput);
    /* check whether fasta is matched */
    if (PathUtil::getFileName(fastaNameInput) != fastaName) {
        fastaNameInput = PathUtil::getAbsPath(fastaNameInput);
        printFastqFileNotMatchInfo(metaRefe);
        LOG_ERROR("initialize reference failed: used fasta %s, should be %s", fastaNameInput.c_str(), fastaName.c_str());
        return false;
    }
    if (fastaLength != fastqFileLenInput) {
        printFastqFileNotMatchInfo(metaRefe);
        LOG_ERROR("initialize reference failed: used fasta file len %ld, should be %ld", fastqFileLenInput, fastaLength);
        return false;
    }

    pRefGene = MemoryUtil::safeNewClass<Reference>(parameter.referenceGenic, parameter.threadNum);
    if (pRefGene == nullptr) {
        return false;
    }
    pRefGene->setNiFile(parameter.niIndexFile);
    /* SAM/BAM need only the reference base sequence, so skip the read-mapping index;
       the type comes from srcFileType, defaulting conservatively to unknown. */
    BlockType srcType = TYPE_UNKNOW;
    {
        const Json::Value* meta = nullptr;
        if (baseFileMeta.getMetaData().isMember("srcFileType")) {
            meta = &baseFileMeta.getMetaData();
        } else if (dynamicFileMeta.getMetaData().isMember("srcFileType")) {
            meta = &dynamicFileMeta.getMetaData();
        }
        if (meta != nullptr) {
            const uint32_t v = (*meta)["srcFileType"].asUInt();
            if (v == (uint32_t)SAM || v == (uint32_t)BAM || v == (uint32_t)FASTQ_GEN2 ||
                v == (uint32_t)FASTQ_GEN3 || v == (uint32_t)BINARY) {
                srcType = (BlockType)v;
            }
        }
    }
    const bool needIndex = !(srcType == SAM || srcType == BAM);
    const bool ok = needIndex ? pRefGene->makeIndex() : pRefGene->makeSquashIndex();
    if (!ok) {
        LOG_ERROR("initialize reference failed");
        MemoryUtil::safeDeleteClass(pRefGene);
        pRefGene = nullptr;
        return false;
    }

    return true;
}

void DecompressEngine::printFastqFileNotMatchInfo(const Json::Value& metaRefe) {
    fprintf(stderr, "need to specify the following FASTA file:\n\n");
    fprintf(stderr, "\t%-12s : %s\n", "File Name", metaRefe["fasta_name"].asString().c_str());
    fprintf(stderr, "\t%-12s : %ld\n", "File Length", metaRefe["fasta_len"].asInt64());
    fprintf(stderr, "\t%-12s : %s\n", "File MD5", metaRefe["fasta_md5"].asString().c_str());
}

bool DecompressEngine::initRefeIndex() {
    std::string indexFileName = parameter.inputFile + ".pbgzi";
    SamIndex::getInstance().loadFromFile(indexFileName);
    return true;
}

int64_t DecompressEngine::readBlocks(BlockReader* blockReader)  {
    if (!parameter.refeGenePos.empty()) {
        return readBlockByPostition(blockReader);
    } else {
        return CodecEngine::readBlocks(blockReader);
    }
}

/*
 * -b conversion stage. The writer thread used to parse every SAM text line and build BAM
 * records itself; with a small -t the conversion budget collapsed to one, so that thread
 * sat at 100% while the coder workers blocked on the output pool. Here the conversion runs
 * on its own thread (mirroring the BAM-read render stage), and the writer only frames the
 * finished records into BGZF.
 */
void DecompressEngine::dispatchOutputBlock(RoughIOBlock* blockPtr) {
    if (bamConvertPool != nullptr) {
        bamConvertPool->push(blockPtr);
    } else {
        outputDataPool->push(blockPtr);
    }
}

int32_t DecompressEngine::startOutputStage() {
    if (!parameter.isDecToBam) {
        return 0;
    }
    bamConvertPool = std::make_unique<BlockingQueueType>();
    /* Keep the in-flight window at least as wide as the worker set, so coder threads do
       not stall waiting for the converter's queue. */
    bamConvertPool->setCapility((parameter.threadNum << 1) + 2);
    bamConvertThreads.emplace_back(&DecompressEngine::bamConvertWorker, this);
    if (parameter.verbose) {
        fprintf(stderr, "[pipeline] SAM->BAM conversion stage up (1 thread)\n");
    }
    return 0;
}

void DecompressEngine::drainOutputStage() {
    if (bamConvertPool == nullptr || bamConvertThreads.empty()) {
        return;
    }
    bamConvertPool->push(nullptr);
    for (auto& t : bamConvertThreads) {
        if (t.joinable()) {
            t.join();
        }
    }
    bamConvertThreads.clear();
}

void DecompressEngine::bamConvertWorker() {
    pthread_setname_np(pthread_self(), "sam2bam");
    /* The per-block line conversion keeps its own small budget; the stage itself is the
       dedicated thread that decouples it from the writer. */
    SamToBamConverter converter(bamConvertThreadNum(parameter.threadNum));
    std::map<int64_t, RoughIOBlock*> pending;
    bool headerDone = false;

    auto forward = [this, &converter](RoughIOBlock* block) {
        if (converter.convert(block) != 0) {
            LOG_ERROR("SAM->BAM conversion failed for block(%ld)", (long)block->getBlockId());
            taskFailed.store(true);
        }
        outputDataPool->push(block);
    };
    auto flushPending = [&]() {
        for (auto& entry : pending) {
            forward(entry.second);
        }
        pending.clear();
    };

    while (true) {
        RoughIOBlock* block = bamConvertPool->get();
        if (block == nullptr) {
            break;
        }

        if (!headerDone) {
            const bool alignment = BlockUtil::isSAMBlock(block->getBlockType()) ||
                                   BlockUtil::isBAMBlock(block->getBlockType());
            /* Block 0 carries the SAM header; a data block can still overtake it in the
               completion queue, so hold such blocks until the header is converted. */
            const bool headerCandidate = block->getBlockId() == 0 ||
                (alignment && block->getDataLen() > 0 && block->getBuffer()[0] == '@');
            if (alignment && !headerCandidate) {
                pending[block->getBlockId()] = block;
                /* Region/partial streams may have no header block; after an in-flight
                   window, synthesize one instead of waiting forever. */
                if (pending.size() >= (size_t)parameter.threadNum + 1) {
                    headerDone = true;
                    auto it = pending.begin();
                    RoughIOBlock* first = it->second;
                    pending.erase(it);
                    forward(first);
                    flushPending();
                }
                continue;
            }
            headerDone = true;
        }

        forward(block);
        if (!pending.empty()) {
            flushPending();
        }
    }
    flushPending();
}

int64_t DecompressEngine::readBlockByPostition(BlockReader* blockReader) {
    if (blockReader == nullptr) {
        return -1;
    }

    // Read the first block to parse the Pbgz file header
    std::unique_ptr<RoughIOBlock> inBlockPtr = std::make_unique<RoughIOBlock>(getBlockSize());
    if (inBlockPtr == nullptr) {
        LOG_ERROR("Get free block failed.");
        return -1;
    }
    std::unique_ptr<RoughIOBlock> outBlockPtr = std::make_unique<RoughIOBlock>(getBlockSize());
    if (outBlockPtr == nullptr) {
        LOG_ERROR("Get free block failed.");
        return -1;
    }

    std::unique_ptr<SamCodecActuator> actuator;
    do {
        inBlockPtr->reset();
        int64_t ret = blockReader->readBlock(inBlockPtr.get(), TYPE_UNKNOW);
        if (ret <= 0) {
            break;
        }
        if (BlockUtil::isAuxiliaryBlock(inBlockPtr->getBlockType())) {
            /* Auxiliary blocks (e.g. the QUAL prior) are not data blocks; hand them to
               the consumer rather than parsing them as a SAM header. */
            PbgzBlockReader* pbgzReader = dynamic_cast<PbgzBlockReader*>(blockReader);
            (void)offerAuxBlock(inBlockPtr.get(),
                                pbgzReader ? pbgzReader->getCurrentFileIndex() : 0);
            continue;
        }
        outBlockPtr->reset();
        actuator = std::make_unique<SamCodecActuator>(inBlockPtr.get(), outBlockPtr.get(), this, pRefGene);
        actuator->initMetaInfo();
        if (actuator->getHeadLineNumber() > 0) {
            actuator->decompressHeader(outBlockPtr.get());
        }
        if (actuator->getHeadLineNumber() == 0 || actuator->getSamLineNumber() > 0) {
            break;
        }
    } while (true);

    uint16_t chrIndex = SamInfo::getInstance().getChrNameIndex(refPosChrName);
    if (chrIndex == 65535) {
        return 0;
    }
    std::set<std::pair<uint32_t, int64_t>> blockList;
    SamIndex::getInstance().getSamBlockByRef(chrIndex, refPosBegin, refPosEnd, blockList);

    FileReader* fileReader = dynamic_cast<FileReader*>(ioReader);
    if (fileReader != nullptr) {
        BlockType fileType = TYPE_UNKNOW;
        for (auto block : blockList) {
            fileReader->seekIO(block.second);
            int64_t ret = readOneBlock(blockReader, fileType);
            if (ret <= 0) {
                break;
            }
        }
    } else {
        PipeReader* pipeReader = dynamic_cast<PipeReader*>(ioReader);
        if (pipeReader != nullptr) {
            BlockType fileType = TYPE_UNKNOW;
            int64_t ret = 0;
            do {
                RoughIOBlock* blockPtr = freeInputPool->get();
                if (blockPtr == nullptr) {
                    LOG_ERROR("Get free block failed.");
                    return -1;
                }
                blockPtr->reset();

                int64_t ret = blockReader->readBlock(blockPtr, fileType);
                if (ret <= 0) {
                    freeInputPool->push(blockPtr);
                    return ret;
                }

                if (BlockUtil::isAuxiliaryBlock(blockPtr->getBlockType())) {
                    PbgzBlockReader* pbgzReader = dynamic_cast<PbgzBlockReader*>(blockReader);
                    (void)offerAuxBlock(blockPtr, pbgzReader ? pbgzReader->getCurrentFileIndex() : 0);
                    freeInputPool->push(blockPtr);
                    continue;
                }

                bool isMatch = false;
                for (auto block : blockList) {
                    if (blockPtr->getBlockId() == block.first) {
                        updateInputStatics(blockPtr);
                        inputDataPool->push(blockPtr);
                        if (ret > 0) {
                            blockCount++;
                        }
                        isMatch = true;
                        break;
                    }
                }
                if (!isMatch) {
                    freeInputPool->push(blockPtr);
                }
            } while (ret > 0 || ret == -2);
        }
    }
    return 0;
}

bool QualPriorConsumer::claim(RoughIOBlock* blockPtr, int64_t packageIndex) {
    if (blockPtr == nullptr || blockPtr->getBlockType() != QUAL_PRIOR) {
        return false;
    }

    {
        /* Already cached (e.g. fetched by a random-access seek): claim without re-decompressing. */
        std::lock_guard<std::mutex> lock(mutex);
        if (byPackage.find(packageIndex) != byPackage.end()) {
            return true;
        }
    }

    coder_json cmeta;
    Json::Value meta;
    cmeta.decoder(blockPtr->getMetaBuffer(), blockPtr->getMetaLen(), meta);
    if (meta["coder"]["magic"].asString() != "bzip2") {
        LOG_ERROR("Unknown qual prior coder: %s", meta["coder"]["magic"].asString().c_str());
        return true;
    }

    const unsigned int srcLen = (unsigned int)meta["srclen"].asUInt64();
    const unsigned int dstLen = (unsigned int)meta["dstlen"].asUInt();
    if (srcLen == 0 || dstLen == 0 || blockPtr->getDataLen() < dstLen) {
        LOG_ERROR("Qual prior block malformed: srclen=%u dstlen=%u datalen=%u",
                  srcLen, dstLen, blockPtr->getDataLen());
        return true;
    }

    std::vector<uint8_t> plain(srcLen);
    unsigned int outLen = srcLen;
    int rc = BZ2_bzBuffToBuffDecompress((char*)plain.data(), &outLen,
                                        (char*)(void*)blockPtr->getBuffer(), dstLen, 0, 0);
    if (rc != BZ_OK || outLen != srcLen) {
        LOG_ERROR("Decompress qual prior failed, rc=%d, outLen=%u, expect=%u", rc, outLen, srcLen);
        return true;
    }

    {
        std::lock_guard<std::mutex> lock(mutex);
        byPackage[packageIndex] = std::make_shared<const std::vector<uint8_t> >(std::move(plain));
    }
    LOG_INFO("Qual prior loaded for package %lld: %u -> %u bytes", (long long)packageIndex, dstLen, srcLen);
    return true;
}

bool DecompressEngine::unpackReference(PbgzBlockReader* blockReader, Json::Value& refeMeta) {
    int64_t refeSquashLen = refeMeta["squash_len"].asInt64();
    std::string fastaName = refeMeta["fasta_name"].asString();
    int64_t maxLen = 16 << 20;
    if (refeMeta.isMember("max_block_len")) {
        maxLen = refeMeta["max_block_len"].asInt64();
    }
    int64_t blocks = refeMeta["blocks"].asInt64();
    std::string md5Packed = refeMeta["md5"].asString();
    int64_t offset = 0;
    if (refeMeta.isMember("offset")) {
        offset = refeMeta["offset"].asInt64();
        LOG_INFO("Reference offset = %lld", offset);
    }

    FileReader* fileReader = nullptr;
    size_t readPos = 0;
    if (offset > 0) {
        fileReader = dynamic_cast<FileReader*>(ioReader);
        if (fileReader == nullptr) {
            return false;
        }
        readPos = fileReader->getCurrentPos();
        /* refe offset is relative to this package's start, not absolute; add the
           package start so concatenated packages still resolve. */
        fileReader->seekIO((int64_t)blockReader->getCurrentFileStart() + offset);
    }
    const int64_t pcnt = parameter.threadNum;
    /* Heap vector because pcnt is runtime-dependent (no VLA). */
    std::vector<RoughIOBlock*> block((size_t)pcnt);
    for (int64_t n = 0; n < pcnt; n++) {
        block[n] = MemoryUtil::safeNewClass<RoughIOBlock>(maxLen);
    }
    
    BlockingQueue<RoughIOBlock*> inputPool(pcnt);
    BlockingQueue<RoughIOBlock*> inputBlock(1);
    pRefGene =  MemoryUtil::safeNewClass<Reference>(fastaName, parameter.threadNum);
    uint8_t* refeSquash = pRefGene->initSquashByStream(refeSquashLen);

    std::vector<std::thread> tpools;
    for (int64_t n = 0; n < blocks; n++) {
        if (n < pcnt) {
            blockReader->readBlock(block[n]);
            inputPool.push(block[n]);
            tpools.push_back(std::thread([&inputPool, &inputBlock, &refeSquash]() {
                for (;;) {
                    RoughIOBlock* currBlock = inputPool.get();
                    if (currBlock == nullptr) {
                        return;
                    }
                    /* Parse this block's meta to get its stream and encoder info */
                    coder_json cmeta;
                    Json::Value refBlockMeta;
                    cmeta.decoder(currBlock->getMetaBuffer(), currBlock->getMetaLen(), refBlockMeta);
                    if (refBlockMeta["coder"]["magic"].asString() == "coder_ppmd") {
                        int64_t offset = refBlockMeta["offset"].asInt64();
                        int32_t srcLen = refBlockMeta["srclen"].asInt();
                        int32_t dstLen = refBlockMeta["dstlen"].asInt();

                        /* Destination offset within the reference squash */
                        if (currBlock->getDataLen() != dstLen) {
                            LOG_ERROR("reference unpack failed in block %ld: encoded length expect %ld, actual is %ld",
                                refBlockMeta["block"].asInt(), dstLen, currBlock->getDataLen());
                            return;
                        }
                        coder_io refeIo(currBlock->getBuffer(), currBlock->getDataLen());
                        refeIo.meta = refBlockMeta["coder"];
                        coder_ppmd cppmd(&refeIo);
                        int32_t len = cppmd.decode(refeSquash + offset, srcLen);
                        if (srcLen != len) {
                            LOG_ERROR("reference unpack failed in block %ld: encoded length expect %ld, actual is %ld",
                                    refBlockMeta["block"].asInt(), srcLen, len);
                            return;
                        } 
                    } else {
                        LOG_ERROR("undefined coder %s", refBlockMeta["coder"]["magic"].asString().c_str());
                        return;
                    }
                    inputBlock.push(currBlock);
                }
            }));
        } else {
            RoughIOBlock* current = inputBlock.get(); 
            blockReader->readBlock(current);
            inputPool.push(current);
        }
    }

    if (fileReader != nullptr && readPos != 0) {
        fileReader->seekIO(readPos);
    }

    uint32_t poolSize = tpools.size();
    for (uint32_t i= 0; i < poolSize; ++i) {
        (void)inputBlock.get();
        inputPool.push(nullptr);
    }

    for (auto &t : tpools) {
        if (t.joinable()) {
            t.join();
        }
    }
    
    for (int64_t n = 0; n < pcnt; n++) {
        MemoryUtil::safeDeleteClass(block[n]);
    }
    return true;
}

Actuator* DecompressEngine::createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
     if (parameter.isMakeIndex) {
        fprintf(stderr, "Binary file will not make index.");
        parameter.isMakeIndex = false;
    }

    Actuator* pActuator = nullptr;
    if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType())) {
        pActuator = MemoryUtil::safeNewClass<FastqDecompressActuator>(inBlockPtr, outBlockPtr, this, pRefGene);
    } else if (BlockUtil::isSAMBlock(inBlockPtr->getBlockType()) || BlockUtil::isBAMBlock(inBlockPtr->getBlockType())) {
        pActuator = MemoryUtil::safeNewClass<SamDecompressActuator>(inBlockPtr, outBlockPtr, this, pRefGene);
    } else if (inBlockPtr->getBlockType() == BINARY) {
        pActuator = MemoryUtil::safeNewClass<BinaryDecompressActuator>(inBlockPtr, outBlockPtr, this);
    }

    if (pActuator == nullptr) {
        /* Don't release blocks here; cleanup is handled uniformly by startWorkTask. */
        LOG_ERROR("Not support block type: %d, blockId=%ld", inBlockPtr->getBlockType(), inBlockPtr->getBlockId());
        return nullptr;
    }

    if (pActuator->initial() != 0) {
        LOG_ERROR("actuator initial failed");
        MemoryUtil::safeDeleteClass(pActuator);
        return nullptr;
    }

    return pActuator;
}
