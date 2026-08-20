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
#include "utils/path_util.h"
#include "coder.h"
#include "coder_json.h"
#include "coder_ppmd.h"
#include "pbgz_types.h"
#include "codec_actuator_adapter.h"

#include <bzlib.h>

BlockReader* DecompressEngine::createBlockReader() {
    PbgzBlockReader* pbgzReader = MemoryUtil::safeNewClass<PbgzBlockReader>(ioReader);
    pbgzReader->init();
    baseFileMeta = pbgzReader->getBaseFileMeta();
    // 读回压缩时写入的块大小上界，供 actuator 预分配输出缓冲（见 PbgzEngine::fileBlockSize）
    if (baseFileMeta.getMetaData().isMember("block_size")) {
        fileBlockSize = baseFileMeta.getMetaData("block_size").asUInt();
    }
    dynamicFileMeta = pbgzReader->getDynamicFileMeta();
    if (dynamicFileMeta.getMetaData().isMember("qual_prior")) {
        (void)unpackQualPrior(pbgzReader, dynamicFileMeta.getMetaData("qual_prior"));
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

            // Parse the p parameter, format like chr1:100-200, convert the part before : to chrID, the part before and after - for the reference gene start position
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
        /* -b：把解压出的 SAM 文本转成标准 BAM 写出 */
        BlockWriter* blockWriter = MemoryUtil::safeNewClass<BamWriter>(ioWriter);
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
    /* -b 的 BamWriter 需要收尾：flush BGZF 残留块并追加 EOF 标记 */
    BamWriter* bamWriter = dynamic_cast<BamWriter*>(blockWriter);
    if (bamWriter != nullptr && 0 != bamWriter->finish()) {
        LOG_ERROR("Bam writer finish failed.");
    }
    MemoryUtil::safeDeleteClass(blockWriter);
}

/*
 * 随机读取用先验：按文件级偏移 seek 过去，再交给块读者解析。
 *
 * 顺序解压靠注册表在路过时认领即可，但随机读会直接跳到任意数据块，它所属包的
 * 先验块可能根本没被读到过，必须能按需取回。偏移记的是块容器头（魔数所在处），
 * 所以这里 seek 完只需要 readBlock，容器的魔数、长度、校验和仍由块读者负责，
 * 不在这里手工解析字节。
 *
 * 与 unpackReference 一样，偏移是相对本包起点的距离而不是绝对位置，cat 之后
 * 必须加上包起点；读完要把读位置还原，否则主循环会从先验块之后接着读。
 */
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
    /* 偏移只用来 seek 取块，登记与查表一律用包序号。 */
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

    /* refe metadata may be in either baseFileMeta or dynamicFileMeta.
       baseFileMeta is written before startWorkPreProc adds refe, so on
       file-to-file compression the refe is only in dynamicFileMeta. */
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
    /*
     * SAM/BAM 解压只取参考碱基序列（getSquash/getStretch 等），不需要 makeIndex 的
     * 读映射哈希表（那只对 FASTQ 的 read 定位有用）。压缩侧把源文件类型写进
     * srcFileType（文件输入在 baseFileMeta，管道输入在 dynamicFileMeta）。
     *
     * 老版本文件、或管道输入压缩的文件可能没有这个字段：getMetaData 会对缺失键造出
     * null，asUInt() 返回 0。这里显式兼容——只有识别为已知类型才采用，否则按未知
     * 保守处理（走 makeIndex）。参考是用户明确下的指令, 用不了就得停下, 悄悄改成
     * 无参考等于解出另一份东西。
     */
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

int64_t DecompressEngine::readBlockByPostition(BlockReader* blockReader) {
    if (blockReader == nullptr) {
        return -1;
    }

    // Read the beginning block of the Pbgz file to complete the file header parsing
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
            /*
             * 头部预读同样要守"辅助块不是数据块"这条规则。先验块现在位于文件首块之前，
             * 拿它去当 SAM 头解析会得到 headLineNumber == 0 而提前跳出循环，
             * 染色体信息就永远注册不上，区域查询直接返回空。
             * 顺手把它交给认领者：这条路径读到的先验与顺序流读到的是同一个块。
             */
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
        /* 缓存已有就直接认领，不重复解压：同一个先验既可能被随机读提前 seek 取回，
           也会在顺序流里被再次路过，两条路径指向的是同一个绝对地址。 */
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
        /*
         * refe 里存的 offset 是相对本压缩包头部的距离，不是拼接文件里的绝对位置。
         *
         * 压缩时输出总是从 0 开始写，所以两者在单包情况下恰好相等，长期以来直接当绝对
         * 位置用也没出过问题。但多个包 cat 拼接之后，第二个包整体后移了前面所有包的
         * 长度，再按绝对位置 seek 就会落进前一个包的数据区，参考解不出来，解压在第一个
         * 包结束处就停住——这正是带参考基因组时 cat 拼接失效的原因。
         *
         * 加上本包起点即可还原。单包时起点为 0，行为与改动前逐字节一致。
         */
        fileReader->seekIO((int64_t)blockReader->getCurrentFileStart() + offset);
    }
    int64_t pcnt = parameter.threadNum;
    RoughIOBlock* block[pcnt];
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
            tpools.push_back(std::thread([&inputPool, &inputBlock, &refeSquash, &n]() {
                for (;;) {
                    RoughIOBlock* currBlock = inputPool.get();
                    if (currBlock == nullptr) {
                        return;
                    }
                    /* First parse out current block's meta information, get corresponding line's stream information and encoder information */
                    coder_json cmeta;
                    Json::Value refBlockMeta;
                    cmeta.decoder(currBlock->getMetaBuffer(), currBlock->getMetaLen(), refBlockMeta);
                    if (refBlockMeta["coder"]["magic"].asString() == "coder_ppmd") {
                        int64_t offset = refBlockMeta["offset"].asInt64();
                        int32_t srcLen = refBlockMeta["srclen"].asInt();
                        int32_t dstLen = refBlockMeta["dstlen"].asInt();

                        /* Point to the position of data to be decompressed in current block of reference */
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
        /* 失败时一律不释放任何块，收尾统一由 startWorkTask 负责，避免两处各归还一次。 */
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
