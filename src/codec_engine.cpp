/*
 * codec_engine.cpp - Cpp file for pbgz project
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

#include "codec_engine.h"
#include "log/logger.h"
#include "coder/coder.h"
#include "utils/memory_util.h"
#include "pbgz_manager.h"
#include "block_wrapper.h"
#include "fastq_actuator.h"
#include "binary_actuator.h"
#include "utils/path_util.h"
#include "coder_ppmd.h"
#include "coder_json.h"
#include "config_manager.h"
#ifdef __SSE4_2__ 
#include "hardware.h"
#endif
#include <set>
#include "sam_actuator.h"

CodecEngine::~CodecEngine() {
    MemoryUtil::safeDeleteClass(pRefGene);
}


void CodecEngine::writeFilePostProc(BlockWriter* blockWriter) {
    PbgzBlockWriter* pbgzWriter =  dynamic_cast<PbgzBlockWriter*>(blockWriter);
    if (pbgzWriter != nullptr) {
        {
            std::lock_guard<std::mutex> lock(dynamicFileMetaMutex);
            if (dynamicFileMeta.getMetaData().empty()) {
                return;
            }
            pbgzWriter->updateHeadExt();
            pbgzWriter->setDynamicFileMeta(dynamicFileMeta);
            pbgzWriter->writeDynamicFileMeta();
        }
        resetReferenceOffset();
    }
    return;
}

void CodecEngine::writeOneBlock(BlockWriter* blockWriter, RoughIOBlock* outBlockPtr) {
    setDataBlockPosition(outBlockPtr->getBlockId());
    return PbgzEngine::writeOneBlock(blockWriter, outBlockPtr);
}

void CodecEngine::writeBlockPreProc(BlockWriter*, RoughIOBlock* outBlockPtr) {
    if (outBlockPtr->getBlockType() == REFERENCE) {
        /*
         * This hook must be called by the writer thread after blockId sorting and
         * before writeBlock. If the offset were captured at outputDataPool enqueue time,
         * the preceding data blocks might not be on disk yet, so refe.offset would vary
         * with worker-thread scheduling and in turn alter the dynamic metadata.
         */
        FileWriter* fileWriter =  dynamic_cast<FileWriter*>(ioWriter);
        if (fileWriter != nullptr) {
            updateReferenceOffset(fileWriter->getCurrentPos());
        }
    }
    return;
}

void CodecEngine::updateReferenceOffset(int64_t offset) {
    std::lock_guard<std::mutex> lock(dynamicFileMetaMutex);
    if (refeOffsetFLag) {
        return;
    }
    LOG_DEBUG("Reference offset is %ld", offset);
    refeOffsetFLag = true;
    dynamicFileMeta.getMetaData("refe")["offset"] = offset;
    return;
}

void CodecEngine::resetReferenceOffset() {
    std::lock_guard<std::mutex> lock(dynamicFileMetaMutex);
    refeOffsetFLag = false;
}
