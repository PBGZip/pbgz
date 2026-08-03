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

int64_t CompressEngine::readOneBlock(BlockReader* blockReader, BlockType& fileType) {
    int64_t ret = CodecEngine::readOneBlock(blockReader, fileType);
    
    // Initialize statistics when file type is determined (first block)
    if (ret > 0 && !statsInitialized && fileType != TYPE_UNKNOW) {
        initStatsBasedOnFileType(fileType);
    }
    
    return ret;
}

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

    runFilePreprocessOnce(inBlockPtr);

    return actuator;
}

void CompressEngine::runFilePreprocessOnce(RoughIOBlock* inBlockPtr) {
    /*
     * runFilePreprocessOnce —— 整个文件生命周期内只跑一次的编码器试压/挑选。
     *
     * 这里的 tryClaim（IDLE→RUNNING 的无锁 CAS）**不是**真正决定"只跑一次"
     * 的地方，它只是最后一道兜底：万一未来有代码路径绕过流水线并发调进
     * 本函数，也能保证 CodecSelector::analyze 不会被重复触发。
     *
     * 真正让首块预处理天然串行的机制在 PbgzEngine::startWorkTask ——
     * 所有 id > 0 的工作线程刚启动就 conditionVar.wait 挂起，只有 id == 0
     * 的线程会走到 createActuator，进而调到本函数，为第一个数据块完成
     * 整套编码器筛选。首块处理完毕后 isNeedNotify 恰好在此刻返回 true，
     * 随后 notify_all 才放行其余线程。也就是说，正常路径上执行到这里时，
     * 其他工作线程都还阻塞在条件变量上；tryClaim 在正常路径永远命中，
     * 抢占逻辑只作为防御性存在。
     *
     * 之所以要把首块预处理死死钉在串行阶段，核心目的是让编码器决策
     * **确定且可复现**：给定同一份输入，被采样的永远是同一段字节，
     * 各字段试压得到的信号也一致，全文件后续所有块选到的 coder 都稳定
     * 不变，压缩产物达到字节级可复现。反之若放任多线程同时抢首块预处理，
     * 各线程看到的样本会随 OS 调度抖动而不同，最终 coder 选择就变成
     * "依赖调度"的非确定行为，回归对比与压测都无从谈起。
     *
     * 预处理内部另有一个**正交**的并发维度："同一份样本、不同 coder 并行
     * 试压"（各 coder 之间相互独立、无副作用），这与流水线"不同数据块、
     * 不同工作线程"的并发维度不是一回事，不要与外层的串并转换混为一谈。
     */
    if (!preprocessInfo.tryClaim()) {
        return;
    }

    if (0 != CodecSelector::analyze(inBlockPtr, preprocessInfo)) {
        LOG_INFO("File preprocessing (codec selection) failed, using default coders");
    } else if (parameter.verbose) {
        static const char* samFieldNames[] = {
            "QNAME", "FLAG", "RNAME", "POS", "MAPQ", "CIGAR",
            "RNEXT", "PNEXT", "TLEN", "SEQ", "QUAL"
        };
        static const char* fqFieldNames[] = {
            "ID", "SEQ", "QUAL", "COMMENT"
        };

        bool isSam = (preprocessInfo.fileType == SAM);
        uint32_t fieldCount = isSam ? SAM_FIELD_COUNT : FQ_FIELD_COUNT;
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
        if (!preprocessInfo.qualPrior().empty()) {
            fprintf(stderr, "  QUAL 先验：训练 %llu 字节，快照 %zu 字节\n",
                    (unsigned long long)preprocessInfo.qualPriorTrainedBytes(),
                    preprocessInfo.qualPrior().size());
        }
        fprintf(stderr, "\n");
    }
    preprocessInfo.markDone();;
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
        if (!pRefGene->makeIndex()) {
            LOG_ERROR("makeIndex failed, reference will not be used");
            MemoryUtil::safeDeleteClass(pRefGene);
            pRefGene = nullptr;
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

    if (BlockUtil::isFastqBlock(inBlockPtr->getBlockType()) || BlockUtil::isSAMBlock(inBlockPtr->getBlockType())) {
        pActuator = actuatorPreProc(pActuator, inBlockPtr, outBlockPtr);
    }

    return pActuator;
}
