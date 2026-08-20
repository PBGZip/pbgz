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
    // Read back the block-size upper bound written during compression, for the actuator to pre-allocate its output buffer (see PbgzEngine::fileBlockSize)
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
        /* -b: convert the decompressed SAM text to standard BAM when writing */
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
    /* The BamWriter under -b needs finalization: flush the residual BGZF block and append the EOF marker */
    BamWriter* bamWriter = dynamic_cast<BamWriter*>(blockWriter);
    if (bamWriter != nullptr && 0 != bamWriter->finish()) {
        LOG_ERROR("Bam writer finish failed.");
    }
    MemoryUtil::safeDeleteClass(blockWriter);
}

/*
 * Prior for random access: seek to the file-level offset, then hand it to the block
 * reader for parsing.
 *
 * Sequential decompression only needs the consumer to claim the prior when it is
 * passed by; random access, however, jumps straight to an arbitrary data block whose
 * package's prior block may never have been read, so the prior must be retrievable on
 * demand. The offset points at the block container header (where the magic number
 * lives), so after seeking, readBlock suffices; the container's magic number, length,
 * and checksum are still handled by the block reader, not parsed manually here.
 *
 * As with unpackReference, the offset is a distance relative to the start of the
 * package, not an absolute position; after concatenating with cat, the package start
 * must be added. The read position must also be restored when done, otherwise the main
 * loop would continue reading from after the prior block.
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
    /* The offset is only used to seek and fetch the block; registration and lookups always use the package index. */
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
     * SAM/BAM decompression only needs the reference base sequence (getSquash/
     * getStretch and so on); it does not need the read-mapping hash table from
     * makeIndex (useful only for locating FASTQ reads). The compression side records
     * the source file type in srcFileType (in baseFileMeta for file input, in
     * dynamicFileMeta for piped input).
     *
     * Old-version files, or files compressed from piped input, may lack this field:
     * getMetaData fabricates null for a missing key and asUInt() returns 0. This is
     * explicitly handled here — the type is adopted only when it is recognized as a
     * known type, otherwise it is treated conservatively as unknown (going through
     * makeIndex). If an explicitly specified reference cannot be loaded, decompression
     * must abort, otherwise it would produce data that does not match expectations.
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
             * Header pre-reading must also obey the rule that "an auxiliary block is
             * not a data block". The prior block now sits before the file's first data
             * block; parsing it as a SAM header would yield headLineNumber == 0 and
             * exit the loop early, so chromosome information would never be registered
             * and region queries would return empty.
             * Hand it to the consumer as well: the prior read on this path is the same
             * block the sequential stream would read.
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
        /* If already cached, claim it directly without decompressing again: the same
           prior may be fetched early by a random-access seek and then encountered again
           in the sequential stream, and both paths point at the same absolute address. */
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
         * The offset stored in refe is a distance relative to the head of this
         * compression package, not an absolute position in a concatenated file.
         *
         * During compression the output always starts writing at 0, so the two happen
         * to be equal in the single-package case, and treating the value as an absolute
         * position has worked for a long time. But once multiple packages are
         * concatenated with cat, the second package is shifted down by the total length
         * of all preceding packages; seeking to the absolute position then lands in the
         * previous package's data region, the reference cannot be unpacked, and
         * decompression stops at the end of the first package — which is exactly why
         * cat concatenation failed when a reference genome was involved.
         *
         * Adding the current package's start restores the address. For a single package
         * the start is 0, so the behavior is byte-for-byte identical to before the
         * change.
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
        /* On failure, never release any block here; cleanup is uniformly handled by startWorkTask to avoid each block being returned twice. */
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
