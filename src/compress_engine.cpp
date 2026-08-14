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

    if (BlockUtil::isSAMBlock(fileType)) {
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

void CompressEngine::releaseBlockReader(BlockReader* &blockReader) {
    MemoryUtil::safeDeleteClass(blockReader);
}

void CompressEngine::releaseBlockWriter(BlockWriter* &blockWriter) {
    MemoryUtil::safeDeleteClass(blockWriter);
}

Actuator* CompressEngine::actuatorPreProc(Actuator* actuator, RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) {
    FastqCompressActuator* fastqActuator = dynamic_cast<FastqCompressActuator*>(actuator);
    if (fastqActuator != nullptr) {
        FastqCodecActuator* fastqCodecActuator = fastqActuator->getCodecActuator();
        if (0 != fastqCodecActuator->preAnalysis()) {
            LOG_INFO("Fastq preAnalysis failed");
            MemoryUtil::safeDeleteClass(fastqActuator);
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
 * 编码器试压/挑选与 QUAL 先验决策，整个文件生命周期内只跑一次。
 *
 * 由读线程在派发第 0 块之前调用（PbgzEngine::readOneBlock）。之所以钉在这个位置：
 *
 * 一是**决策必须先于任何块的压缩**。放在入队之前，这条先后关系由数据流位置保证，
 *   不再依赖"其余工作线程还挂在条件变量上"这种时序巧合。
 * 二是**决策结果必须确定可复现**。读线程是唯一的数据源头，被采样的永远是第 0 块的
 *   同一段字节，各字段试压信号一致，全文件后续所有块选到的 coder 稳定不变，
 *   压缩产物达到字节级可复现。此前若首块的 preAnalysis 失败，本函数会改由某个并行
 *   工作线程在它自己那一块上触发，采样内容随调度漂移——那才是真正的隐患。
 * 三是**读线程手上有文件级信息**（输入总长、可否 seek），先验是否划算这类按全文件
 *   数据量做的判断才有立足点。
 *
 * 注意：这里只做**决策**。QUAL 先验的训练和发布被推迟到跨块预训练完成
 * （finalizePretrain）——首块通常只有 ~9 MB QUAL，跨块累积到 45 MB 才能逼近
 * 先验收益的收敛点。发布期间 coder 线程在 workStartBarrier 门闩上等待。
 *
 * tryClaim（IDLE→RUNNING 的无锁 CAS）现在纯属兜底：调用方只有读线程一个，且只在
 * 第 0 块调，重入不可能发生。
 *
 * 预处理内部另有一个**正交**的并发维度："同一份样本、不同 coder 并行试压"
 * （各 coder 相互独立、无副作用），与流水线"不同数据块、不同工作线程"不是一回事。
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
     * 文件总长只有读线程手上有，也只有在这里才谈得上"按全文件摊销划不划算"。
     * 管道输入拿不到长度，返回 0 表示不可知，由判决器自己决定怎么处理。
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

        bool isSam = (preprocessInfo.fileType == SAM);
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
                 * 候选是谁由 FieldCodecSelection 自己声明，这里照着念即可。
                 * 早先是两个写死的槽位加上"质量值行换标签"的补丁，候选超过两个之后
                 * 那种写法就不成立了。吞吐按 样本字节 / 试压耗时 折算成 MB/s。
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
                fprintf(stderr, "  %-6s -> %-14s  %u -> %u (%.2f%%)  [%s] %u轮\n",
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
     * 先验值得训练：进入跨块预训练。块 0 由随后的 pretrainBlockProc 累积 QUAL，
     * 读完目标块数（=并发度）或 45 MB 上限后由 finalizePretrain 统一训练、发布，
     * 届时才通知 coder 线程开始压缩。此处先不 emit、不 markDone。
     */
    if (preprocessInfo.wantQualPrior()) {
        pretrainTarget = (parameter.threadNum > 0) ? parameter.threadNum : 1;
        pretraining = true;
        pretrainBlockCount = 0;
        return;
    }

    /* 无先验：保持原行为——发射（空先验时为 no-op）、标记完成、放行 coder 线程。 */
    emitQualPrior();
    preprocessInfo.markDone();
    signalPriorSettled();

    if (parameter.verbose) {
        fprintf(stderr, "[流水线] 编码器决策完成，开始并行压缩（%u 线程）\n", parameter.threadNum);
    }
}

/*
 * 逐块累积 QUAL 先验训练样本。块 0 的决策已在 fileDecisionProc 完成，这里对所有
 * 预训练块（含块 0）追加 QUAL，累积到目标块数或 45 MB 上限后训练并发布。
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
 * 训练累积的 QUAL、把先验作为辅助块发布落盘、置 DONE，并通知等待中的 coder 线程。
 * 调用点：预训练块数够/45 MB 满（读循环中），或文件读完仍不足时（readLoopPostProc）。
 * 数据块要引用先验，先验就必须先于一切数据块被压缩完成；发布在本函数内同步完成。
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
        qualParams = qualSel->fcv2Params;   /* 与压缩端同档位训练，先验才用得上 */
    }
    std::vector<uint8_t> snapshot =
        CodecSelector::trainQualPriorModel(qualPriorAccum, qualParams, &trainedBytes);
    qualPriorAccum = {};   /* 释放训练样本内存 */

    if (!snapshot.empty()) {
        preprocessInfo.setQualPrior(std::move(snapshot), trainedBytes);
        if (parameter.verbose) {
            fprintf(stderr, "  QUAL 先验：跨块训练 %llu 字节，快照 %zu 字节\n",
                    (unsigned long long)trainedBytes,
                    preprocessInfo.qualPrior().size());
        }
    }

    emitQualPrior();
    preprocessInfo.markDone();
    signalPriorSettled();

    if (parameter.verbose) {
        fprintf(stderr, "[流水线] 先验发布完成，开始并行压缩（%u 线程）\n", parameter.threadNum);
    }
}

/* 读循环结束兜底：预训练未完成（文件不足目标块数）则收尾；无先验/空文件也放行。 */
void CompressEngine::readLoopPostProc() {
    if (pretraining) {
        finalizePretrain();
    }
    signalPriorSettled();
}

/* coder 线程启动门闩：等读线程发布先验（或判定无先验）的通知。 */
void CompressEngine::workStartBarrier() {
    std::unique_lock<std::mutex> lk(priorMutex);
    priorCv.wait(lk, [this] { return priorSettled; });
}

void CompressEngine::emitQualPrior() {
    const std::vector<uint8_t>& prior = preprocessInfo.qualPrior();
    if (prior.empty()) {
        return;
    }

    /* 管道输出无法 seek，解压侧取不到先验，写了也是净损失。 */
    FileWriter* fileWriter = dynamic_cast<FileWriter*>(ioWriter);
    if (fileWriter == nullptr) {
        return;
    }

    /* bzip2 的最坏膨胀是 1% + 600 字节，官方接口要求调用方按此备足。 */
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
     * 不给辅助块分配数据块 id：它位置寻址，写线程见到即写。
     * 早先这里做的是 blockId2Write++，那是把写线程的游标当作发号器用，
     * 在数据块开写之前发射就会直接顶掉第 0 块的位置。
     */
    outBlock->setBlockId(-1);

    /* 落盘与取址都发生在写线程内，返回后 outBlock 已被归还，不可再触碰。 */
    const int64_t offset = emitSyncAuxBlock(outBlock);
    if (offset < 0) {
        LOG_ERROR("Emit qual prior block failed");
        return;
    }
    /* 先发布内容再发布地址：取用方以地址原子量的 acquire 作为两者的可见性关口。 */
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
        /* 参考是用户明确下的指令, 用不了就得停下, 悄悄改成无参考等于压出另一份东西。 */
        if (!pRefGene->makeIndex()) {
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
                if (currBlockLen >= outBlock->getRemain()) {
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
  * 只负责把参考基因组的描述填进基础文件元信息，不产生任何 IO。
  *
  * 这段内容原先和下面 startWorkPreProc 里的 packReference 写在一起，而
  * startWorkPreProc 是在写线程启动之后才执行的。写线程一起来就会把 baseFileMeta
  * 落盘，于是"主线程填 refe"和"写线程写 baseFileMeta"变成了两个线程对同一份
  * JSON 的赛跑：写线程赢，文件头里就没有 refe；主线程赢，就有。
  *
  * 后果是同一份输入压两次得到的字节不同（实测两个稳定变体相差 128 字节，正好是
  * refe 这个 JSON 成员压缩后的长度），后续所有内容整体偏移。两个结果都能正确解压，
  * 因为解压侧对 refe 缺失做了兼容（见 decompress_engine.cpp 里关于 refe 可能落在
  * baseFileMeta 或 dynamicFileMeta 的说明），所以这个问题一直没有暴露出来。
  *
  * 现在把纯元信息的构造提前到写线程启动之前，竞争消失；真正需要写线程就绪才能做的
  * packReference 留在 startWorkPreProc 不动。
  */
 int32_t CompressEngine::prepareFileMeta() {
    /*
     * 块大小是压缩时输入分块的上界，写进文件头让解压侧按它预分配 block_size*2
     * 的输出缓冲——确定上界，不是估算，不受 fieldcount/读长影响。这是堵所有
     * 字段越界的主防线；coder_io 的 putc 检查与 decode 错误返回链作兜底。
     */
    baseFileMeta.setMetaData("block_size", getBlockSize());

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
  * 参考基因组打包在文件头的场景下，把参考块写往下游。
  *
  * 这一步必须留在写线程启动之后，因为 packReference 会把参考块推进输出队列。
  * 元信息本身已经在 prepareFileMeta 里填好，这里不再触碰 baseFileMeta——
  * 一旦在这里改它，就会重新引入上面描述的那个竞争。
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
    } else if (BlockUtil::isSAMBlock(inBlockPtr->getBlockType())) {
        pActuator = MemoryUtil::safeNewClass<SamCompressActuator>(inBlockPtr, outBlockPtr, this, pRefGene);
    } else if (inBlockPtr->getBlockType() == BINARY) {
        pActuator = MemoryUtil::safeNewClass<BinaryCompressActuator>(inBlockPtr, outBlockPtr, this);
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

    if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType()) || BlockUtil::isSAMBlock(inBlockPtr->getBlockType())) {
        pActuator = actuatorPreProc(pActuator, inBlockPtr, outBlockPtr);
    }

    return pActuator;
}
