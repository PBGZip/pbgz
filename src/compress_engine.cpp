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
#include "utils/path_util.h"
#include "coder_ppmd.h"
#include "coder_json.h"
#include "codec_actuator_adapter.h"
#include "codec_selector.h"
#include "field_coder_config.h"

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
    BlockReader* blockReader = BlockFactory::createBlockReader(ioReader, parameter.compressLevel);
    if (blockReader == nullptr) {
        LOG_ERROR("Create block reader failed.");
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
    /* Pre-analysis is done on the coder side: state (idSplitSymbols/contentPos/qualFreqTable etc.) stays in the actuator for use during compression */
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
 * Codec trial compression/selection and QUAL prior decision, run only once over the
 * entire file lifetime.
 *
 * Called by the reader thread before dispatching block 0 (PbgzEngine::readOneBlock).
 * It is pinned to this position for three reasons:
 *
 * First, the decision must precede the compression of any block. Placing it before
 * the block is enqueued guarantees this ordering through the data-flow position,
 * rather than relying on a timing coincidence such as "the remaining worker threads
 * are still waiting on a condition variable".
 * Second, the decision must be deterministic and reproducible. The reader thread is
 * the sole source of data, so the sampled bytes are always the same segment of block 0;
 * the trial-compression signal for every field is identical, the coder selected for
 * all subsequent blocks is stable across the file, and the compressed output is
 * reproducible at the byte level. Previously, if the preAnalysis of the first block
 * failed, this function would instead be triggered by some parallel worker thread on
 * its own block, so the sampled content drifted with scheduling — that was the real
 * hazard.
 * Third, only the reader thread holds file-level information (total input length,
 * whether the input is seekable), which grounds whole-file judgments such as whether
 * a prior is worthwhile.
 *
 * Note: only the decision is made here. QUAL prior training and publication are
 * deferred until cross-block pretraining completes (finalizePretrain) — the first
 * block usually carries only ~9 MB of QUAL, and accumulating across blocks up to
 * 45 MB is needed to approach the convergence point of the prior's benefit. While the
 * prior is being published, coder threads wait at the workStartBarrier latch.
 *
 * tryClaim (an IDLE→RUNNING lock-free CAS) is now purely a safety net: the only
 * caller is the reader thread, and it is invoked only on block 0, so re-entrancy
 * cannot occur.
 *
 * Preprocessing has another orthogonal concurrency dimension internally:
 * "trial-compressing the same sample with different coders in parallel" (each coder
 * is independent and side-effect free), which is not the same as the pipeline's
 * "different data blocks on different worker threads".
 */
void CompressEngine::fileDecisionProc(RoughIOBlock* inBlockPtr) {
    CodecEngine::fileDecisionProc(inBlockPtr);

    if (!statsInitialized) {
        initStatsBasedOnFileType(inBlockPtr->getBlockType());
    }
    if (!preprocessInfo.tryClaim()) {
        return;
    }

    /*
     * Only the reader thread knows the total file length, and only here does it make
     * sense to ask "is it worth amortizing across the whole file". Piped input has no
     * length available; returning 0 means unknown, and the decision logic decides how
     * to handle it.
     */
    FileReader* fileReader = dynamic_cast<FileReader*>(ioReader);
    const uint64_t inputTotalBytes = (fileReader != nullptr && fileReader->getFileSize() > 0)
                                     ? (uint64_t)fileReader->getFileSize() : 0;

    if (0 != CodecSelector::analyze(inBlockPtr, inputTotalBytes, preprocessInfo)) {
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
                /*
                 * The candidate coders are declared by FieldCodecSelection itself; here
                 * we simply echo them. Previously this was two hardcoded slots plus a
                 * patch that relabeled the quality-value line; that approach broke once
                 * there were more than two candidates. Throughput is derived as sample
                 * bytes / trial-compression time, in MB/s.
                 */
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

    /*
     * For piped input the source type cannot be detected in prepareFileMeta, so
     * baseFileMeta has no srcFileType. By this point the type of the first data block
     * is known, so it is recorded in the dynamic metadata; the decompressor then does
     * not need to build the read-mapping index.
     */
    if (parameter.inputFile == STDIN) {
        const BlockType t = preprocessInfo.fileType;
        if (t == SAM || t == BAM || t == FASTQ_GEN2 || t == FASTQ_GEN3 || t == BINARY) {
            std::lock_guard<std::mutex> lock(dynamicFileMetaMutex);
            dynamicFileMeta.setMetaData("srcFileType", (Json::Value::UInt)t);
        }
    }

    /*
     * The prior is worth training: enter cross-block pretraining. Block 0 accumulates
     * QUAL via the subsequent pretrainBlockProc; once the target block count
     * (= concurrency) or the 45 MB cap is reached, finalizePretrain trains and
     * publishes it, and only then notifies the coder threads to start compressing.
     * Nothing is emitted or marked done here.
     */
    if (preprocessInfo.wantQualPrior()) {
        pretrainTarget = (parameter.threadNum > 0) ? parameter.threadNum : 1;
        pretraining = true;
        pretrainBlockCount = 0;
        return;
    }

    /* No prior: keep the original behavior — emit (a no-op when the prior is empty), mark done, and release the coder threads. */
    emitQualPrior();
    preprocessInfo.markDone();
    signalPriorSettled();

    if (parameter.verbose) {
        fprintf(stderr, "[pipeline] coder decision complete, starting parallel compression (%u threads)\n", parameter.threadNum);
    }
}

/*
 * Accumulate QUAL prior training samples block by block. The block-0 decision is
 * already made in fileDecisionProc; here QUAL is appended for every pretraining block
 * (including block 0), and once the target block count or the 45 MB cap is reached the
 * prior is trained and published.
 */
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

/*
 * Train on the accumulated QUAL, publish the prior to disk as an auxiliary block,
 * set DONE, and notify the waiting coder threads.
 * Call sites: once the pretraining block count is reached / 45 MB is full (in the read
 * loop), or when the file ends before enough blocks were read (readLoopPostProc).
 * Because data blocks reference the prior, the prior must be fully compressed before
 * any data block; publication completes synchronously within this function.
 */
void CompressEngine::finalizePretrain() {
    if (!pretraining) {
        return;
    }
    pretraining = false;

    uint64_t trainedBytes = 0;
    QualFcv2Params qualParams;
    const FieldCodecSelection* qualSel =
        (preprocessInfo.fields.size() > (size_t)SAM_QUAL)
            ? &preprocessInfo.fields[SAM_QUAL] : nullptr;
    if (qualSel != nullptr && qualSel->selectedCoder == CoderType::FCV2) {
        qualParams = qualSel->fcv2Params;   /* train with the same parameter tier as the compression side, otherwise the prior is unusable */
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

/* Read-loop fallback at end: finalize if pretraining did not complete (fewer blocks than target in the file); files without a prior / empty files are also released. */
void CompressEngine::readLoopPostProc() {
    if (pretraining) {
        finalizePretrain();
    }
    signalPriorSettled();
}

/* coder thread startup latch: waits for the reader thread's notification that the prior is published (or decided absent). */
void CompressEngine::workStartBarrier() {
    std::unique_lock<std::mutex> lk(priorMutex);
    priorCv.wait(lk, [this] { return priorSettled; });
}

void CompressEngine::emitQualPrior() {
    const std::vector<uint8_t>& prior = preprocessInfo.qualPrior();
    if (prior.empty()) {
        return;
    }

    /* Piped output cannot be seeked, so the decompressor could never retrieve the prior; writing it would be a net loss. */
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    if (fileWriter == nullptr) {
        return;
    }

    /* bzip2's worst-case expansion is 1% + 600 bytes; the official API requires callers to size the buffer accordingly. */
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
    /*
     * The auxiliary block gets no data-block id: it is addressed by position, and the
     * writer thread writes it as soon as it is seen. Previously this did
     * blockId2Write++, using the writer thread's cursor as an id generator; emitting
     * before any data block had been written would directly overwrite the position of
     * block 0.
     */
    outBlock->setBlockId(-1);

    /* Both the write-to-disk and the address retrieval happen on the writer thread; after this returns, outBlock has been returned to the pool and must not be touched. */
    const int64_t offset = emitSyncAuxBlock(outBlock);
    if (offset < 0) {
        LOG_ERROR("Emit qual prior block failed");
        return;
    }
    /* Publish the content before the address: consumers use the acquire load of the address atomic as the visibility barrier for both. */
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
        /*
         * SAM/BAM carry their own alignment positions (RNAME/POS); compressing SEQ only
         * needs the reference base sequence, so the read-mapping hash table from
         * makeIndex (useful only for locating FASTQ reads) is unnecessary. Piped input
         * has no type available, so it falls back to the conservative index-building
         * path. If an explicitly specified reference cannot be loaded, compression must
         * abort, otherwise it would produce data that does not match expectations.
         */
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
    BlockingQueueType* freeOutput = freeOutputPool.get();
    BlockingQueueType* outputPool = outputDataPool.get();
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
                if ((size_t)currBlockLen >= outBlock->getRemain()) {
                    LOG_ERROR("reference block pack failed");
                    break;
                }
                memcpy(outBlock->getCurrent(), refeIo.data, refeIo.data_len);
                outBlock->setDataLen(refeIo.data_len);
                memcpy(outBlock->getCurrent(), metaString.c_str(), metaString.length());
                outBlock->setMetaLen(metaString.length());
                outBlock->setBlockType(REFERENCE);
                outBlock->setBlockId(count + refe2do.first.first);
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

 /*
  * This only fills the reference genome description into the base file metadata; it
  * performs no I/O.
  *
  * This content used to be written together with packReference in startWorkPreProc
  * below, but startWorkPreProc runs only after the writer thread has started. As soon
  * as the writer thread starts it writes baseFileMeta to disk, so "the main thread
  * filling in refe" and "the writer thread writing baseFileMeta" became two threads
  * racing over the same JSON: if the writer thread won, the file header had no refe;
  * if the main thread won, it did.
  *
  * The consequence was that compressing the same input twice produced different bytes
  * (two stable variants were measured to differ by 128 bytes, exactly the compressed
  * length of the refe JSON member), shifting everything downstream. Both results
  * decompress correctly, because the decompressor tolerates a missing refe (see the
  * note in decompress_engine.cpp about refe possibly living in baseFileMeta or
  * dynamicFileMeta), which is why the problem never surfaced.
  *
  * Now the construction of this pure metadata is moved before the writer thread starts,
  * eliminating the race; packReference, which genuinely needs the writer thread to be
  * ready, stays in startWorkPreProc.
  */
 int32_t CompressEngine::prepareFileMeta() {
    /*
     * The block size is the upper bound for splitting the input during compression; it
     * is written into the file header so the decompressor can pre-allocate a
     * block_size*2 output buffer. This is a definite upper bound, not an estimate,
     * unaffected by field count / read length. It is the primary defense against every
     * field overflowing its buffer; coder_io's putc checks and the decode error-return
     * chain act as fallbacks.
     */
    baseFileMeta.setMetaData("block_size", getBlockSize());

    /*
     * The source file type is written into the file header so the decompressor can
     * decide whether the reference only needs to load its squash (SAM/BAM) or must
     * build the read-mapping index (FASTQ). Piped input has no type available and it
     * is not written (the decompressor then treats the type as unknown,
     * conservatively).
     */
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
    refeMeta["fasta_md5"] = pRefGene->getFastaChecksum();
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

 /*
  * In the scenario where the reference genome is packed in the file header, push the
  * reference blocks downstream.
  *
  * This step must stay after the writer thread starts, because packReference pushes
  * reference blocks into the output queue. The metadata itself is already filled in by
  * prepareFileMeta; baseFileMeta must not be touched here — modifying it here would
  * reintroduce the race described above.
  */
 int32_t CompressEngine::startWorkPreProc() {
    if (pRefGene == nullptr || !isPackRefeInHeader()) {
        return 0;
    }

    int64_t maxRefLen = 0;
    int64_t totolEncLen = 0;
    (void)packReference(maxRefLen, totolEncLen, false);

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
    } else if (BlockUtil::isSAMBlock(inBlockPtr->getBlockType()) || BlockUtil::isBAMBlock(inBlockPtr->getBlockType())) {
        pActuator = MemoryUtil::safeNewClass<SamCompressActuator>(inBlockPtr, outBlockPtr, this, pRefGene);
    } else if (inBlockPtr->getBlockType() == BINARY) {
        pActuator = MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr, this);
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

    if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType()) || BlockUtil::isSAMBlock(inBlockPtr->getBlockType())
       || BlockUtil::isBAMBlock(inBlockPtr->getBlockType())) {
        pActuator = actuatorPreProc(pActuator, inBlockPtr, outBlockPtr);
    }

    return pActuator;
}
