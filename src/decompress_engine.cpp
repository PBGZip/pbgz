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

BlockReader* DecompressEngine::createBlockReader() {
    PbgzBlockReader* pbgzReader = MemoryUtil::safeNewClass<PbgzBlockReader>(ioReader);
    pbgzReader->init();
    baseFileMeta = pbgzReader->getBaseFileMeta();
    dynamicFileMeta = pbgzReader->getDynamicFileMeta();
    if (!initRefGene(pbgzReader)) {
        LOG_INFO("Init reference for decompress failed");
        MemoryUtil::safeDeleteClass(pbgzReader);
        return nullptr;
    }

    if (parameter.inputFile != STDIN && !parameter.refeGenePos.empty()) {
        if (!initRefeIndex()) {
            LOG_INFO("Init reference index for decompress failed");
        }
    }

    return pbgzReader;
}

BlockWriter* DecompressEngine::createBlockWriter() {
    BlockWriter* blockWriter = MemoryUtil::safeNewClass<BlockWriter>(ioWriter);
    if (blockWriter == nullptr) {
        LOG_ERROR("Failed to create block writer.");
        return nullptr;
    }
    blockWriter->init();
    return blockWriter;
}

void DecompressEngine::releaseBlockReader(BlockReader* &blockReader) {
    MemoryUtil::safeDeleteClass(blockReader);
}

void DecompressEngine::releaseBlockWriter(BlockWriter* &BlockWriter) {
    MemoryUtil::safeDeleteClass(BlockWriter);
}

bool DecompressEngine::initRefGene(PbgzBlockReader* blockReader) {
    if (blockReader == nullptr) {
        return false;
    }
    if (!baseFileMeta.getMetaData("refe").isObject()) {
        LOG_DEBUG("No reference");
        return true;
    }
    Json::Value& metaRefe = baseFileMeta.getMetaData("refe");
    if (dynamicFileMeta.getMetaData().isMember("refe")) {
        LOG_DEBUG("Read meta data from dynamic file meta");
        metaRefe = dynamicFileMeta.getMetaData("refe");
    }

    std::string fastaName = metaRefe["fasta_name"].asString();
    int64_t fastaLength = metaRefe["fasta_len"].asInt64();
    std::string niName = metaRefe["ni_name"].asString();

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

    Reference refeCheck(parameter.referenceGenic, parameter.threadNum);
    std::string niNameInput;
    refeCheck.getNiFileFromReference(niNameInput);
    niNameInput = PathUtil::getAbsPath(niNameInput);
    if (PathUtil::getFileName(niNameInput) != niName)  {
        LOG_ERROR("initialize reference failed: used ni file %s, should be %s", niNameInput.c_str(), niName.c_str());
        return false;
    }
    /* matched, do make index */
    pRefGene = MemoryUtil::safeNewClass<Reference>(parameter.referenceGenic, parameter.threadNum);
    if (pRefGene == nullptr) {
        return false;
    }
    if (!pRefGene->initSquashByNiFile()) {
        LOG_ERROR("initialize reference failed");
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
    if (pbgzIndex.initIndex() != 0) {
        return false;
    }
    return true;
}

void DecompressEngine::readBlocks(BlockReader* blockReader)  {
    if (!parameter.refeGenePos.empty() && pRefGene) {
        return readBlockByPostition(blockReader);
    } else {
        return CodecEngine::readBlocks(blockReader);
    }
}

void DecompressEngine::readBlockByPostition(BlockReader* blockReader) {
    if (blockReader == nullptr) {
        return;
    }
    return CodecEngine::readBlocks(blockReader);
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
        fileReader->seekIO(offset);
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
        pActuator = MemoryUtil::safeNewClass<FastqDecompressActuator>(inBlockPtr, outBlockPtr, pRefGene);
    } else if (BlockUtil::isSAMBlock(inBlockPtr->getBlockType())) {
        pActuator = MemoryUtil::safeNewClass<SamDecompressActuator>(inBlockPtr, outBlockPtr, pRefGene);
    } else if (inBlockPtr->getBlockType() == BINARY) {
        pActuator = MemoryUtil::safeNewClass<BinaryDecompressActuator>(inBlockPtr, outBlockPtr);
    }
    
    if (pActuator == nullptr) {
        LOG_ERROR("Not support block type: %d, blockId=%d", inBlockPtr->getBlockType(), inBlockPtr->getBlockId());
        freeInputPool->push(inBlockPtr);
        outBlockPtr->reset();
        outBlockPtr->setBlockId(inBlockPtr->getBlockId());
        // When an error occurs, push a block with length 0 but correct ID, the write thread ignores blocks with length 0 to prevent thread waiting
        outputDataPool->push(outBlockPtr);
        return nullptr;
    }

    if (pActuator->initial() != 0) {
        LOG_ERROR("actuator initial failed");
        MemoryUtil::safeDeleteClass(pActuator);
        return nullptr;
    }

    return pActuator;
}
