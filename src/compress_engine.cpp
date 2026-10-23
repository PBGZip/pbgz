/*
 * compress_engine.cpp - Source file for pbgz project
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


#include "compress_engine.h"
#include "profile_stats.h"
#include "utils/md5_util.h"
#include "utils/path_util.h"
#include "coder_ppmd.h"
#include "coder_json.h"
#include "codec_actuator_adapter.h"
#include "codec_selector.h"
#include "field_coder_config.h"
#include "pbgz_index.h"   /* SamIndex / BlockPosition, used by the index paths below */
#include "sam_info.h"

#include <bzlib.h>

int32_t CompressEngine::init() {
    if (0 != CodecEngine::init()) {
        return -1;
    }

    indexBlockQueue->setCapility(2);
    freeIndexBlockQueue->setCapility(2);
    for (uint32_t i = 0; i < freeIndexBlockQueue->getCapility(); ++i) {
        RoughIOBlock* ptr = MemoryUtil::safeNewClass<RoughIOBlock>();
        if (ptr == nullptr) {
            LOG_ERROR("Failed to create RoughIOBlock for freeIndexBlockQueue");
            return -1;
        }
        freeIndexBlockQueue->push(ptr);
    }

    return 0;
}

void CompressEngine::initStatsBasedOnFileType(BlockType fileType) {
    if (statsInitialized) {
        return;
    }

    if (BlockUtil::isSAMBlock(fileType) || BlockUtil::isBAMBlock(fileType)) {
        stats = std::make_unique<SamStat>();
        if (stats->init() != 0) {
            LOG_ERROR("Failed to initialize SamStat");
        } else {
            LOG_INFO("Initialized SamStat for SAM file compression");
        }
    } else if (BlockUtil::isFastqBlock(fileType)) {
        stats = std::make_unique<FastqStat>();
        if (stats->init() != 0) {
            LOG_ERROR("Failed to initialize FastqStat");
        } else {
            LOG_INFO("Initialized FastqStat for FASTQ file compression");
        }
    } else {
        LOG_DEBUG("File type %d does not require stats initialization", fileType);
    }

    statsInitialized = true;
}

CompressEngine::~CompressEngine() {
    while(!freeIndexBlockQueue->empty()) {
        RoughIOBlock* ptr = freeIndexBlockQueue->get();
        MemoryUtil::safeDeleteClass(ptr);
    }

    while(!indexBlockQueue->empty()) {
        RoughIOBlock* ptr = indexBlockQueue->get();
        MemoryUtil::safeDeleteClass(ptr);
    }
}

BlockReader* CompressEngine::createBlockReader() {
    const bool bamStruct = (parameter.mode == PBGZ_MODE_FAST);
    /* pbgzAsOpaque: read a pbgz file as opaque binary instead of decoding it. bgzfThreads
       sizes the parallel BAM inflate stage (one per six coding threads, see
       bamInflateThreadNum). */
    BlockReader* blockReader = BlockFactory::createBlockReader(
        ioReader, parameter.compressLevel, true, bamStruct, true, bamInflateThreadNum(parameter.threadNum));
    if (blockReader == nullptr) {
        LOG_ERROR("Create block reader failed.");
        return nullptr;
    }
    /* Archive mode: stream raw BAM records after the leading text blocks so BAM->SAM
       rendering and MD5 run on the render thread pool instead of the reader. */
    if (!bamStruct) {
        if (BamGzBlockReader* bam = dynamic_cast<BamGzBlockReader*>(blockReader)) {
            const int64_t textLimit = (parameter.threadNum > 0) ? (int64_t)parameter.threadNum : 1;
            bam->setRawStreamMode(true, textLimit);
        } else if (BamBlockReader* bam = dynamic_cast<BamBlockReader*>(blockReader)) {
            const int64_t textLimit = (parameter.threadNum > 0) ? (int64_t)parameter.threadNum : 1;
            bam->setRawStreamMode(true, textLimit);
        }
    }
    return blockReader;
}

BlockWriter* CompressEngine::createBlockWriter() {
    PbgzBlockWriter* blockWriter = BlockFactory::createPbgzBlockWriter(ioWriter);
    if (blockWriter == nullptr) {
        LOG_ERROR("pbgzWriter is NULL.");
        return nullptr;
    }
    blockWriter->setBaseFileMeta(baseFileMeta);
    blockWriter->writeBaseFileMeta();
    return blockWriter;
}

void CompressEngine::releaseBlockReader(BlockReader* &blockReader) {
    MemoryUtil::safeDeleteClass(blockReader);
}

void CompressEngine::releaseBlockWriter(BlockWriter* &blockWriter) {
    MemoryUtil::safeDeleteClass(blockWriter);
}

Actuator* CompressEngine::actuatorPreProc(Actuator* actuator, RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
    /* Pre-analysis runs on the coder side; state stays in the actuator for compression. */
    FastqCompressActuator* fastqActuator = dynamic_cast<FastqCompressActuator*>(actuator);
    if (fastqActuator != nullptr) {
        FastqCodecActuator* fastqCodecActuator = fastqActuator->getCodecActuator();
        if (0 != fastqCodecActuator->preAnalysis()) {
            LOG_INFO("Fastq preAnalysis failed");
            MemoryUtil::safeDeleteClass(fastqActuator);
            inBlockPtr->setBlockType(BINARY);
            BinaryCompressActuator* binaryActuator = MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr, this);
            if (binaryActuator == nullptr) {
                return nullptr;
            }
            binaryActuator->initial();
            return binaryActuator;
        }
    }

    SamCompressActuator* samActuator = dynamic_cast<SamCompressActuator*>(actuator);
    if (samActuator != nullptr) {
        SamCodecActuator* samCodecActuator = samActuator->getCodecActuator();
        if (0 != samCodecActuator->preAnalysis()) {
            LOG_INFO("sam preAnalysis failed");
            MemoryUtil::safeDeleteClass(samActuator);
            inBlockPtr->setBlockType(BINARY);
            BinaryCompressActuator* binaryActuator = MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr, this);
            if (binaryActuator == nullptr) {
                return nullptr;
            }
            binaryActuator->initial();
            return binaryActuator;
        }
    }

    return actuator;
}

/*
 * Codec trial / QUAL prior decision, run once per file by the reader thread before
 * block 0 is dispatched. Prior training is deferred to finalizePretrain.
 */
/*
 * Per-block input-only work done before the block is queued: block MD5 and SAM header
 * chromosome pre-parse.
 */
void CompressEngine::preDispatchBlock(RoughIOBlock* blockPtr) {
    if (blockPtr == nullptr) {
        return;
    }
    const BlockType type = blockPtr->getBlockType();
    if (!BlockUtil::isSAMBlock(type) && !BlockUtil::isFastqBlock(type) && type != BINARY) {
        return;
    }
    if (BlockUtil::isSAMBlock(type)) {
        /* Build the chromosome table on the reader thread before the block is queued,
           so coder threads find it ready (see SamUtil::prefetchHeaderChromosomes). */
        PBGZ_PROF_SCOPE(pbgzprof::READ_HEADER);
        SamUtil::prefetchHeaderChromosomes(blockPtr->getBuffer(), blockPtr->getDataLen(),
                                           blockPtr->getNpos());
    }
    if (blockPtr->getBamColumns() != nullptr) {
        /* Structured BAM: hashing is deferred to the compressing worker. */
        return;
    }
    if (blockPtr->isRawBam()) {
        /* Archive-mode raw-record block: render and MD5 are done by the render pool. */
        return;
    }
    PBGZ_PROF_SCOPE(pbgzprof::READ_MD5);
    std::string md5;
    calcMd5sum(md5, blockPtr->getBuffer(), (uint32_t)blockPtr->getDataLen());
    blockPtr->setMd5(md5);
}

/* Archive-mode BAM render stage: BAM->SAM conversion and per-block MD5 run on a
   dedicated thread pool, off the reader and the compressing workers. */
void CompressEngine::engineReaderReady(BlockReader* blockReader) {
    if (renderStageEnabled) {
        return;
    }
    if (parameter.mode == PBGZ_MODE_FAST) {
        return;
    }
    const bool isBam =
        dynamic_cast<BamGzBlockReader*>(blockReader) != nullptr ||
        dynamic_cast<BamBlockReader*>(blockReader) != nullptr;
    if (!isBam) {
        return;
    }
    renderThreadNum = bamConvertThreadNum(parameter.threadNum);
    renderStageEnabled = true;
    renderInputPool = std::make_unique<BlockingQueueType>();
    renderInputPool->setCapility(parameter.threadNum + 1);
    for (uint32_t i = 0; i < renderThreadNum; ++i) {
        renderThreads.emplace_back(std::thread(&CompressEngine::renderWorker, this));
    }
    if (parameter.verbose) {
        fprintf(stderr, "[pipeline] SAM render stage up (%u conversion thread(s))\n",
                renderThreadNum);
    }
}

void CompressEngine::dispatchInputBlock(RoughIOBlock* blockPtr) {
    if (renderStageEnabled && blockPtr->isRawBam() && renderInputPool != nullptr) {
        renderInputPool->push(blockPtr);
        return;
    }
    inputDataPool->push(blockPtr);
}

void CompressEngine::renderStageDrain() {
    if (!renderStageEnabled) {
        return;
    }
    renderStageEnabled = false;
    for (uint32_t i = 0; i < renderThreadNum; ++i) {
        renderInputPool->push(nullptr);
    }
    for (auto& t : renderThreads) {
        if (t.joinable()) {
            t.join();
        }
    }
    renderThreads.clear();
    renderInputPool.reset();
}

void CompressEngine::renderWorker() {
    pthread_setname_np(pthread_self(), "samrender");
    while (true) {
        RoughIOBlock* inBlockPtr = renderInputPool->get();
        if (inBlockPtr == nullptr) {
            break;
        }
        PBGZ_PROF_SCOPE(pbgzprof::READ_RENDER);
        if (0 != bamRenderRawBlockToText(inBlockPtr)) {
            LOG_ERROR("Render raw BAM block(%ld) to SAM text failed.",
                      (long)inBlockPtr->getBlockId());
            /* Block is left raw on failure, so the worker's isRawBam() fallback
               retries the render or downgrades it to BINARY. */
        }
        inputDataPool->push(inBlockPtr);
    }
}

void CompressEngine::fileDecisionProc(RoughIOBlock* inBlockPtr) {
    CodecEngine::fileDecisionProc(inBlockPtr);

    if (!statsInitialized) {
        initStatsBasedOnFileType(inBlockPtr->getBlockType());
    }
    if (!preprocessInfo.tryClaim()) {
        return;
    }

    /* Total input length in the units the samples were taken in (uncompressed); 0 means
       unknown (piped input). Only the reader thread knows it. */
    uint64_t inputTotalBytes = 0;
    FileReader* fileReader = dynamic_cast<FileReader*>(ioReader);
    if (fileReader != nullptr && fileReader->getFileSize() > 0) {
        inputTotalBytes = (uint64_t)fileReader->getFileSize();
    } else if (parameter.inputFile != STDIN) {
        const int64_t uncompressed = PathUtil::getGzUncompressedSize(parameter.inputFile);
        if (uncompressed > 0) {
            inputTotalBytes = (uint64_t)uncompressed;
            LOG_INFO("Input length from gzip trailer: %lld bytes uncompressed.", (long long)uncompressed);
        }
    }

    /* -m selects the per-field candidate lists; CoderFactory additionally gates them
       by profile/level (see CodecSelector::analyze). */
    const int32_t analyzeRet = CodecSelector::analyze(inBlockPtr, inputTotalBytes, preprocessInfo,
                                                      parameter.compressLevel, parameter.mode,
                                                      getReference());

    /* Reference is passed down so SEQ is trialled on the payload it will be written as
       (see CodecSelector::selectSeqReferenceCoder). */

    if (0 != analyzeRet) {
        LOG_INFO("File preprocessing (codec selection) failed, using default coders");
    } else if (parameter.verbose) {
        static const char* samFieldNames[] = {
            "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR",
            "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL", "OPTION"
        };
        static const char* fqFieldNames[] = {
            "ID", "SEQ", "QUAL", "COMMENT"
        };

        bool isSam = BlockUtil::isSAMBlock(preprocessInfo.fileType) || BlockUtil::isBAMBlock(preprocessInfo.fileType);
        uint32_t fieldCount = isSam ? SAM_FIELD_COUNT_SELECT : FQ_FIELD_COUNT;
        const char* const* names = isSam ? samFieldNames : fqFieldNames;

        fprintf(stderr, "\n[preprocess] codec selection"
                " (scanned %u raw bytes, field samples total %u):\n",
                preprocessInfo.scannedBytes, preprocessInfo.sampleBytes);
        for (uint32_t i = 0; i < fieldCount && i < preprocessInfo.fields.size(); ++i) {
            const auto& sel = preprocessInfo.fields[i];
            const char* status = (sel.status == FieldStatus::SELECTED) ? "selected" :
                                 (sel.status == FieldStatus::SKIPPED) ? "skipped" : "failed";
            if (sel.status == FieldStatus::SELECTED) {
                /* Echo the candidates recorded by FieldCodecSelection; throughput is
                   sample bytes / trial time, in MB/s. */
                std::string trials;
                for (uint32_t t = 0; t < sel.trialCount; ++t) {
                    char one[128];
                    const double mbps = sel.trialUs[t]
                        ? (double)sel.decidedLen / sel.trialUs[t] : 0.0;
                    snprintf(one, sizeof(one), "%s%s %.2f%% @%.0fMB/s",
                             (t == 0) ? "" : " | ",
                             coderTypeToMagic(sel.trialCoder[t]),
                             sel.decidedLen ? 100.0 * sel.trialLen[t] / sel.decidedLen : 0.0,
                             mbps);
                    trials += one;
                }
                fprintf(stderr, "  %-6s -> %-14s  %u -> %u (%.2f%%)  [%s] %u rounds\n",
                        names[i], coderTypeToMagic(sel.selectedCoder),
                        sel.decidedLen, sel.bestCompLen,
                        sel.decidedLen ? 100.0 * sel.bestCompLen / sel.decidedLen : 0.0,
                        trials.c_str(),
                        sel.rounds);
            } else {
                fprintf(stderr, "  %-6s -> %-14s  (sample %u bytes)\n",
                        names[i], status, sel.sampleLen);
            }
        }
        fprintf(stderr, "\n");
    }

    /* Piped input has no source type in baseFileMeta; record the first block's type in
       the dynamic metadata instead. */
    if (parameter.inputFile == STDIN) {
        const BlockType t = preprocessInfo.fileType;
        if (t == SAM || t == BAM || t == FASTQ_GEN2 || t == FASTQ_GEN3 || t == BINARY) {
            std::lock_guard<std::mutex> lock(dynamicFileMetaMutex);
            dynamicFileMeta.setMetaData("srcFileType", (Json::Value::UInt)t);
        }
    }

    /* POS delta prior: produced by analyze, written hex-encoded into the file meta
       (only ~512 bytes, so no auxiliary block needed). */
    {
        const std::vector<uint8_t>& posPrior = preprocessInfo.posPrior();
        if (!posPrior.empty()) {
            static const char hexd[] = "0123456789abcdef";
            std::string hex;
            hex.reserve(posPrior.size() * 2);
            for (uint8_t b : posPrior) {
                hex.push_back(hexd[b >> 4]);
                hex.push_back(hexd[b & 0xf]);
            }
            std::lock_guard<std::mutex> lock(dynamicFileMetaMutex);
            dynamicFileMeta.setMetaData("pos_prior", hex);
            if (parameter.verbose) {
                fprintf(stderr, "  POS prior: %zu bytes written to file meta\n",
                        posPrior.size());
            }
        }
    }

    /* Prior is worth training: enter cross-block pretraining; finalizePretrain later
       trains and publishes it and releases the coder threads. */
    if (preprocessInfo.wantQualPrior()) {
        pretrainTarget = (parameter.threadNum > 0) ? parameter.threadNum : 1;
        pretraining = true;
        pretrainBlockCount = 0;
        return;
    }

    /* No prior: emit (no-op when empty), mark done, release the coder threads. */
    emitQualPrior();
    preprocessInfo.markDone();
    signalPriorSettled();

    if (parameter.verbose) {
        fprintf(stderr, "[pipeline] coder decision complete, starting parallel compression (%u threads)\n", parameter.threadNum);
    }
}

/* Append a pretraining block's QUAL to the prior accumulator; train and publish once
   the target block count or the 45 MB cap is reached. */
void CompressEngine::pretrainBlockProc(RoughIOBlock* blockPtr) {
    if (!pretraining) {
        return;
    }
    ++pretrainBlockCount;
    CodecSelector::accumulateQualPrior(blockPtr, qualPriorAccum);
    if (pretrainBlockCount >= pretrainTarget || qualPriorAccum.full()) {
        finalizePretrain();
    }
}

/* Train on the accumulated QUAL, publish the prior, and notify waiting coder threads.
   Runs from the read loop or readLoopPostProc; publication completes synchronously. */
void CompressEngine::finalizePretrain() {
    PBGZ_PROF_SCOPE(pbgzprof::READ_PRETRAIN_TRAIN);
    if (!pretraining) {
        return;
    }
    pretraining = false;

    uint64_t trainedBytes = 0;
    QualFcv2Params qualParams;
    /* The fcv2 tier comes from SAM_QUAL or FQ_QUAL depending on the format. */
    const uint32_t qualFieldIdx = BlockUtil::isFastqBlock(preprocessInfo.fileType)
                                      ? (uint32_t)FQ_QUAL : (uint32_t)SAM_QUAL;
    const FieldCodecSelection* qualSel =
        (preprocessInfo.fields.size() > (size_t)qualFieldIdx)
            ? &preprocessInfo.fields[qualFieldIdx] : nullptr;
    if (qualSel != nullptr && qualSel->selectedCoder == CoderType::FCV2) {
        qualParams = qualSel->fcv2Params;   /* same tier as the compression side, or the prior is unusable */
    }
    std::vector<uint8_t> snapshot =
        CodecSelector::trainQualPriorModel(qualPriorAccum, qualParams, &trainedBytes);
    qualPriorAccum = {};   /* release training-sample memory */

    if (!snapshot.empty()) {
        preprocessInfo.setQualPrior(std::move(snapshot), trainedBytes);
        if (parameter.verbose) {
            fprintf(stderr, "  QUAL prior: %llu bytes trained across blocks, snapshot %zu bytes\n",
                    (unsigned long long)trainedBytes,
                    preprocessInfo.qualPrior().size());
        }
    }

    emitQualPrior();
    preprocessInfo.markDone();
    signalPriorSettled();

    if (parameter.verbose) {
        fprintf(stderr, "[pipeline] prior published, starting parallel compression (%u threads)\n", parameter.threadNum);
    }
}

/* End-of-read fallback: finalize incomplete pretraining, then release coder threads. */
void CompressEngine::readLoopPostProc() {
    if (pretraining) {
        finalizePretrain();
    }
    signalPriorSettled();
}

/* Coder thread startup latch: waits until the prior is published or decided absent. */
void CompressEngine::workStartBarrier() {
    std::unique_lock<std::mutex> lk(priorMutex);
    priorCv.wait(lk, [this] { return priorSettled; });
}

void CompressEngine::emitQualPrior() {
    const std::vector<uint8_t>& prior = preprocessInfo.qualPrior();
    if (prior.empty()) {
        return;
    }

    /* Piped output is not seekable, so the prior could never be retrieved. */
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    if (fileWriter == nullptr) {
        return;
    }

    /* bzip2 worst-case expansion is 1% + 600 bytes. */
    unsigned int dstLen = (unsigned int)(prior.size() + prior.size() / 100 + 600);
    std::vector<uint8_t> comp(dstLen);
    int rc = BZ2_bzBuffToBuffCompress((char*)comp.data(), &dstLen,
                                      (char*)(void*)prior.data(),
                                      (unsigned int)prior.size(),
                                      9, 0, 30);
    if (rc != BZ_OK) {
        LOG_ERROR("Compress qual prior failed, bzip2 rc = %d", rc);
        return;
    }

    RoughIOBlock* outBlock = freeOutputPool->get();
    if (outBlock == nullptr) {
        LOG_ERROR("Get free block for qual prior failed");
        return;
    }
    outBlock->reset();

    Json::Value meta;
    meta["srclen"] = (Json::Value::UInt64)prior.size();
    meta["dstlen"] = (Json::Value::UInt)dstLen;
    meta["coder"]["magic"] = "bzip2";
    coder_json cmeta;
    std::string metaString;
    cmeta.encoder(meta, metaString);

    if ((uint32_t)(dstLen + metaString.length()) >= outBlock->getRemain()) {
        LOG_ERROR("Qual prior block too large: %u + %zu", dstLen, metaString.length());
        freeOutputPool->push(outBlock);
        return;
    }

    memcpy(outBlock->getCurrent(), comp.data(), dstLen);
    outBlock->setDataLen(dstLen);
    memcpy(outBlock->getCurrent(), metaString.c_str(), metaString.length());
    outBlock->setMetaLen(metaString.length());
    outBlock->setBlockType(QUAL_PRIOR);
    /* Auxiliary block gets no data-block id; it is addressed by position. */
    outBlock->setBlockId(-1);

    /* After this returns, outBlock is back in the pool and must not be touched. */
    const int64_t offset = emitSyncAuxBlock(outBlock);
    if (offset < 0) {
        LOG_ERROR("Emit qual prior block failed");
        return;
    }
    /* Publish content before address; consumers acquire-load the address as the barrier. */
    qualPriorBlob = std::make_shared<const std::vector<uint8_t> >(prior);
    qualPriorOffset.store(offset, std::memory_order_release);

    Json::Value priorMeta;
    priorMeta["offset"] = (Json::Value::Int64)offset;
    priorMeta["srclen"] = (Json::Value::UInt64)prior.size();
    priorMeta["dstlen"] = (Json::Value::UInt)dstLen;
    priorMeta["trained_bytes"] = (Json::Value::UInt64)preprocessInfo.qualPriorTrainedBytes();
    {
        std::lock_guard<std::mutex> lock(dynamicFileMetaMutex);
        dynamicFileMeta.setMetaData("qual_prior", priorMeta);
    }

    LOG_INFO("Qual prior packed: %zu -> %u bytes (%.2f%%) at offset %ld",
             prior.size(), dstLen,
             (double)dstLen * 100.0 / (double)prior.size(), (long)offset);
}

void CompressEngine::writeFilePostProc(BlockWriter* blockWriter) {
    CodecEngine::writeFilePostProc(blockWriter);
    // make index
    if (parameter.isMakeIndex) {
        std::string indexFileName = parameter.outputFile + ".pbgzi";
        if (PathUtil::fileExists(indexFileName)) {
            PathUtil::removeFile(indexFileName);
        }
        SamIndex::getInstance().updateFileOffsetsFromBlockPosition();
        SamIndex::getInstance().dumpToFile(indexFileName);
    }
}


int32_t CompressEngine::startEnginePreProc() {
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
        pRefGene->setNiFile(parameter.niIndexFile);
        /* SAM/BAM need only the reference base sequence, so skip the read-mapping
           index; piped input falls back to the conservative makeIndex path. */
        const bool needIndex =
            (parameter.inputFile == STDIN) ||
            !BlockUtil::isSAMBlock(BlockUtil::detectInputFileType(parameter.inputFile));
        const bool ok = needIndex ? pRefGene->makeIndex() : pRefGene->makeSquashIndex();
        if (!ok) {
            LOG_ERROR("Build reference index from %s failed.", parameter.referenceGenic.c_str());
            MemoryUtil::safeDeleteClass(pRefGene);
            pRefGene = nullptr;
            return false;
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
int64_t CompressEngine::packReference(int64_t &maxBlockLen, int64_t &totalEncLen, bool isSanitizeRef,
                                      bool syncAux) {
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
    BlockingQueueType* freeOutput = freeOutputPool.get();
    BlockingQueueType* outputPool = outputDataPool.get();
    uint32_t count = blockCount;

    // Start threads
    for (uint32_t idx = 0; idx < pcnt; ++idx) {
        tpools.push_back(std::thread([&inputPool, &m, &freeOutput, &outputPool, &refe, &each, &maxBlockLen, &totalEncLen, &isSanitizeRef, &count, syncAux]() {
            uint8_t* output = MemoryUtil::safeAlloc<uint8_t>(each);
            while(true) {
                RefeInfo refe2do = inputPool.get();
                int64_t plen = refe2do.second.second;
                if (plen == 0) {
                    break;
                }

                RoughIOBlock* outBlock = freeOutput->get();
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

                /* Write compressed stream of current block's meta */
                outBlock->reset();
                int64_t currBlockLen = metaString.length() + refeIo.data_len;
                /* A 16 MB reference chunk can compress past the pool block's initial size
                   (which follows -l), so grow it instead of failing the pack. */
                if (outBlock->ensureCapacity((size_t)currBlockLen + 1) != 0 ||
                    (size_t)currBlockLen >= outBlock->getRemain()) {
                    LOG_ERROR("reference block pack failed");
                    break;
                }
                memcpy(outBlock->getCurrent(), refeIo.data, refeIo.data_len);
                outBlock->setDataLen(refeIo.data_len);
                memcpy(outBlock->getCurrent(), metaString.c_str(), metaString.length());
                outBlock->setMetaLen(metaString.length());
                outBlock->setBlockType(REFERENCE);
                outBlock->setBlockId(count + refe2do.first.first);
                /* Header packing (stdout) must reach the writer before any data block; sync-aux
                   blocks bypass block-id ordering, whose ids would otherwise collide with data. */
                outBlock->setSyncAux(syncAux);
                outputPool->push(outBlock);

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

int32_t CompressEngine::startEnginePostProc() {
    bool isPackRefeInTail = pRefGene && !parameter.isUnpackRef && parameter.outputFile != STDOUT;
    if (isPackRefeInTail) {
        int64_t maxRefLen = 0;
        int64_t totolEncLen = 0;
        Json::Value refeMeta;
        refeMeta["squash_len"] = (Json::Value::Int64)(pRefGene->getSquashLength());
        refeMeta["fasta_name"] = PathUtil::getFileName(pRefGene->getFastaFileName());
        refeMeta["fasta_len"] = (Json::Value::Int64)(PathUtil::getFileSize(pRefGene->getFastaFileName()));
        refeMeta["fasta_md5"] = pRefGene->getFastaChecksum();
        refeMeta["ni_name"] = "";
        refeMeta["offset"]= 0;
        {
            std::lock_guard<std::mutex> lock(dynamicFileMetaMutex);
            dynamicFileMeta.setMetaData("refe", refeMeta);
        }
        FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
        if (fileWriter != nullptr) {
            int64_t refBlockCount = packReference(maxRefLen, totolEncLen);
            std::lock_guard<std::mutex> lock(dynamicFileMetaMutex);
            dynamicFileMeta.getMetaData("refe")["max_block_len"] = maxRefLen;
            dynamicFileMeta.getMetaData("refe")["blocks"] = refBlockCount;
        }
    }

    return CodecEngine::engineStartAfterProc();
}

void CompressEngine::setDataBlockPosition(uint32_t blockId) {
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    if (fileWriter == nullptr) {
        return;
    }
    BlockPosition::getInstance().setBlockPosition(blockId, fileWriter->getCurrentPos());
    return;
 }

 /* Fill the reference genome description into base file metadata (no I/O); must run
    before the writer thread starts to avoid a metadata race. */
 int32_t CompressEngine::prepareFileMeta() {
    /* Block size: definite upper bound on a compressed block, written to the header so
       the decompressor can pre-allocate its output buffer. */
    baseFileMeta.setMetaData("block_size", getBlockSize());

    /* Source file type: lets the decompressor decide between loading the reference
       squash (SAM/BAM) and building the read-mapping index (FASTQ). */
    if (parameter.inputFile != STDIN) {
        const BlockType srcType = BlockUtil::detectInputFileType(parameter.inputFile);
        if (srcType == BINARY || srcType == FASTQ_GEN2 || srcType == FASTQ_GEN3 ||
            srcType == SAM || srcType == BAM) {
            baseFileMeta.setMetaData("srcFileType", (Json::Value::UInt)srcType);
        }
    }

    if (pRefGene == nullptr) {
        return 0;
    }

    Json::Value refeMeta;
    refeMeta["squash_len"] = (Json::Value::Int64)(pRefGene->getSquashLength());
    refeMeta["fasta_name"] = PathUtil::getFileName(pRefGene->getFastaFileName());
    refeMeta["fasta_len"] = (Json::Value::Int64)(PathUtil::getFileSize(pRefGene->getFastaFileName()));
    /* The FASTA checksum pass runs in the background against the squash build. When the archive's
       refe meta is written to the end-of-file dynamic meta (file output, not --refunpack) it is
       filled in there after compression, so the load need not block on it; otherwise (stdout, or
       --refunpack) the header must carry it and this waits for the pass. */
    const bool deferFastaMd5 = !parameter.isUnpackRef && parameter.outputFile != STDOUT;
    refeMeta["fasta_md5"] = deferFastaMd5 ? "" : pRefGene->getFastaChecksum();
    refeMeta["ni_name"] = "";
    refeMeta["max_block_len"] = 16 << 20;
    if (isPackRefeInHeader()) {
        refeMeta["blocks"] = calcPackRefeBlockSize();
    } else {
        refeMeta["blocks"] = 0;
    }
    baseFileMeta.setMetaData("refe", refeMeta);

    return 0;
 }

 /* When the reference is packed in the header, push its blocks downstream. Must stay
    after the writer thread starts (packReference fills the output queue). */
 int32_t CompressEngine::startWorkPreProc() {
    if (pRefGene == nullptr || !isPackRefeInHeader()) {
        return 0;
    }

    int64_t maxRefLen = 0;
    int64_t totolEncLen = 0;
    (void)packReference(maxRefLen, totolEncLen, false, /* syncAux */ true);

    return 0;
 }

Actuator* CompressEngine::createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
    Actuator* pActuator = nullptr;
    if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType())) {
        if (parameter.isMakeIndex) {
            fprintf(stderr, "Fastq file will not make index.");
            parameter.isMakeIndex = false;
        }
        pActuator = MemoryUtil::safeNewClass<FastqCompressActuator>(inBlockPtr, outBlockPtr, this, pRefGene);
    } else if (inBlockPtr->getBamColumns() != nullptr) {
        /* Structured BAM block (PBGZ_BAM_STRUCT=1): columns instead of SAM text. */
        pActuator = MemoryUtil::safeNewClass<BamCompressActuator>(inBlockPtr, outBlockPtr, this, pRefGene);
    } else if (BlockUtil::isSAMBlock(inBlockPtr->getBlockType()) || BlockUtil::isBAMBlock(inBlockPtr->getBlockType())) {
        pActuator = MemoryUtil::safeNewClass<SamCompressActuator>(inBlockPtr, outBlockPtr, this, pRefGene);
    } else if (inBlockPtr->getBlockType() == BINARY) {
        pActuator = MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr, this);
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

    if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType()) || BlockUtil::isSAMBlock(inBlockPtr->getBlockType())
       || BlockUtil::isBAMBlock(inBlockPtr->getBlockType())) {
        pActuator = actuatorPreProc(pActuator, inBlockPtr, outBlockPtr);
    }

    return pActuator;
}
