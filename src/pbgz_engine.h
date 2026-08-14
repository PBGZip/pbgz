/*
 * pbgz_engine.h - Head file for pbgz project
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

#pragma once

#include <thread>
#include <list>
#include <atomic>

#include "pbgz_types.h"
#include "blocking_queue.h"
#include "io_block.h"
#include "io_wrapper.h"
#include "block_wrapper.h"
#include "actuator.h"
#include "utils/timer.h"
#include "reference.h"
#include "aux_block_consumer.h"

struct PreprocessInfo;

using BlockingQueueType = BlockingQueue<RoughIOBlock*>;

class PbgzEngine {
public:
    PbgzEngine(const PbgzParameter&  para);

    virtual ~PbgzEngine();

    virtual int32_t init();

    virtual int32_t start();

protected:
    virtual void printHeadInfo() { };

    virtual void printTailInfo(Timer&) { }

    virtual uint32_t getBlockSize();

    virtual int32_t startEnginePreProc() { return 0; }

    /*
     * 准备基础文件元信息，必须在 startWriteTask 之前调用。
     *
     * 为什么要单独开这个钩子：写线程一启动就会执行 createBlockWriter，里面立刻
     * 调用 writeBaseFileMeta 把 baseFileMeta 落盘。如果此时主线程还在往
     * baseFileMeta 里塞内容，两边就在同一份 JSON 上赛跑——谁先跑完决定了写进
     * 文件头的元信息里有没有那部分内容。
     *
     * 这个竞争的实际后果是压缩结果不可复现：同一份输入压两次，文件头元信息可能
     * 相差一整个 JSON 成员，后续所有内容随之整体偏移，两次输出的字节完全不同。
     * 两个结果都能正确解压（解压侧对成员缺失做了兼容），所以问题长期没有暴露，
     * 但它让"同一输入必得同一输出"这条基本性质不成立，也让任何字节级回归对比
     * 都无法进行。
     *
     * 因此约定：凡是要写进基础文件元信息的内容，都必须在本钩子里填好；
     * 需要写线程就绪之后才能做的事（例如往下游写参考块）留在 startWorkPreProc。
     */
    virtual int32_t prepareFileMeta() { return 0; }

    virtual int32_t startWriteTask();

    virtual int32_t startWorkPreProc() { return 0; }

    virtual int32_t startWorkTask();

    virtual int32_t startReadPreProc() { return 0; }

    virtual int32_t startReadTask();

    virtual int32_t startEnginePostProc() { return 0; }

    virtual BlockReader* createBlockReader() = 0;

    virtual BlockWriter* createBlockWriter() = 0;

    virtual void releaseBlockReader(BlockReader* &blockReader) = 0;

    virtual void releaseBlockWriter(BlockWriter* &blockWriter) = 0;

    virtual Actuator* createActuator(RoughIOBlock* inBlockPtr, RoughIOBlock* outBlockPtr) = 0;

    virtual int32_t actuatorProc(Actuator* actuator, RoughIOBlock*, RoughIOBlock*);

    virtual int64_t readOneBlock(BlockReader* blockReader, BlockType& fileType);

    virtual int64_t readBlocks(BlockReader* blockReader);

    /*     * 读到一个数据块之后、入队之前回调，供引擎做按块累积等预处理相关工作。
     *
     * 基类默认空操作。压缩引擎用它跨块累积 QUAL 先验训练样本（首块决策 + 后续块
     * 追加），累积到目标块数或上限后训练并发布先验。
     */
    virtual void pretrainBlockProc(RoughIOBlock* /*blockPtr*/) {}

    /*     * 读循环结束（EOF 或读错误）后调用一次，供引擎收尾未完成的预处理。
     * 基类默认空操作。
     */
    virtual void readLoopPostProc() {}

    /*     * 工作线程开始拉块压缩前的同步点。基类默认空操作。
     *
     * 压缩引擎用它实现"先验发布前 coder 线程不得开始"：读线程把跨块训练出的先验
     * 发布（落盘辅助块 + 置位标志）后通知，工作线程才被放行。其他引擎无此需求。
     */
    virtual void workStartBarrier() {}

    /*     * 读到一个块之后、送进 inputDataPool 之前的处置，由 readOneBlock 统一询问。
     *
     * 之所以做成一个返回值而不是让各引擎重写整个 readOneBlock：读循环的骨架
     * （取空闲块 / 读 / 首块决策 / 入队 / 计数）对所有引擎完全相同，各引擎真正
     * 不同的只有"这个块归不归我"这一步。此前三个引擎各抄了一份骨架，任何加在
     * 骨架上的改动都会漏掉两份。
     */
    enum class BlockIntake {
        DISPATCH,   /* 正常数据块，入队交给工作线程 */
        SKIP,       /* 不属于本引擎的数据流，已就地消费，读循环继续 */
        ABORT       /* 输入不可用，读循环终止 */
    };

    virtual BlockIntake intakeBlock(BlockReader* /*blockReader*/, RoughIOBlock* /*blockPtr*/) {
        return BlockIntake::DISPATCH;
    }

    /* 读一个块之前的时机，用于捕捉"这个块在源文件里的起始位置"这类只在读前有效的量。 */
    virtual void readBlockPreProc(BlockReader* /*blockReader*/) { }

    /*     * 文件级决策：读线程在把第一个数据块送进 inputDataPool 之前调一次。
     *
     * 位置是这个钩子的全部意义。编码器选型、先验训练这类决策对整个文件只做一次，
     * 且必须先于任何块的处理。放在"入队之前"，这条先后关系就由数据流位置本身保证：
     * 工作线程只能从队列里取块，取不到就无块可压。队列的入队/出队自带 release/acquire，
     * 决策结果的可见性一并解决。因此不需要工作线程互相等待，也不需要额外的同步标志。
     *
     * 一并得到的性质：决策永远发生在读线程、永远基于第 0 块，样本不再随调度漂移；
     * 同步辅助块（如 QUAL 先验）发射时写线程尚未收到任何数据块，"辅助块物理上先于
     * 全部数据块"从时序上的巧合变成位置上的事实。
     */
    virtual void fileDecisionProc(RoughIOBlock* /*firstBlock*/) { }

    virtual void updateInputStatics(RoughIOBlock*) { }

    virtual void updateOutputStatics(RoughIOBlock*) { }

    virtual void writeBlockPreProc(BlockWriter*, RoughIOBlock*) { }

    virtual void writeBlockPostProc() { }

    virtual void writeFilePreProc() { }

    virtual void writeFilePostProc(BlockWriter*) { }

    virtual void writeOneBlock(BlockWriter* blockWriter, RoughIOBlock* outblockPtr);

    /* 写线程侧：落盘一个同步辅助块并把它的容器头偏移回交给发射者。 */
    void completeSyncAuxBlock(BlockWriter* blockWriter, RoughIOBlock* block);

    /*
     * 路过一个辅助块时逐个询问认领者。返回值仅用于日志：没人认领的辅助块被安静跳过，
     * 这是老版本读到新格式时需要的前向兼容行为，不构成错误。
     */
    bool offerAuxBlock(RoughIOBlock* blockPtr, int64_t packageIndex);

    /*
     * 同步发射一个辅助块：推给写线程，阻塞等它落盘，返回该块**容器头**的绝对文件偏移
     * （即块读管理器能从该处解析出块 json 的位置，不是块内负载的起点）。
     *
     * 辅助块是位置寻址的：它不占用数据块的 blockId，也不参与写线程的顺序重排，
     * 写线程见到即写。之所以必须由写线程亲自执行写动作，是因为只有它持有 BlockWriter
     * 并独占文件写指针——调用方在别的线程上取 getCurrentPos() 拿到的都是竞态值。
     *
     * 仅允许在 fileDecisionProc 内调用：那时读线程还没派发过任何数据块，写线程手上
     * 必然空闲，全局至多一个发射在途，故下面用单个槽位承接结果，无需按块建映射。
     * 返回后块已被写线程归还到 freeOutputPool，调用方不得再触碰它。
     */
    int64_t emitSyncAuxBlock(RoughIOBlock* block);

    void registerAuxConsumer(AuxBlockConsumer* consumer) {
        if (consumer != nullptr) {
            auxConsumers.push_back(consumer);
        }
    }

public:
    std::unique_ptr<BlockingQueueType> freeInputPool;   // Free queue, file reading tasks get blocks from here to read data
    std::unique_ptr<BlockingQueueType> inputDataPool;   // Data reading tasks write completed data to this queue, compression/decompression tasks get data from this queue
    std::unique_ptr<BlockingQueueType> freeOutputPool;  // Compression/decompression tasks get blocks from this queue to write processed data, output tasks write free blocks to this queue after processing
    std::unique_ptr<BlockingQueueType> outputDataPool;  // After compression/decompression is completed, write to this queue, output tasks get data from this queue
    PbgzParameter parameter;
    IOReader* ioReader;
    IOWriter* ioWriter;
    std::list<RoughIOBlock*> outputSortedCache;
    const PbgzParameter& getParameter() const { return parameter; }
    virtual Reference* getReference() { return nullptr; }

    /*
     * 压缩时写入文件头的块大小上界。解压侧 createBlockReader 从 baseFileMeta 读回，
     * actuator 据此 ensureCapacity(block_size*2) 预分配输出缓冲——堵所有字段越界的
     * 主防线。0 表示尚未读到（老文件没写或 createBlockReader 未跑），actuator 落回
     * 默认 getBlockSize()。coder_io 的 putc 检查与 decode 错误返回链作兜底。
     */
    uint32_t fileBlockSize = 0;
    uint32_t getFileBlockSize() const { return fileBlockSize; }

    /* File preprocessing result (codec pre-selection). Only the compression
       engine populates this; other engines return nullptr. */
    virtual const PreprocessInfo* getPreprocessInfo() { return nullptr; }

    /*
     * 取 QUAL 先验模型快照。两侧来源不同但接口一致：
     * 压缩侧返回首块预处理训练出的那份；解压侧按辅助块绝对地址取回已认领的那份。
     * 返回空表示本次没有可用先验，调用方必须回退到固定初值模型。
     *
     * 解压侧不做 seek 回读：先验块现在物理上位于全部数据块之前，
     * 顺序流一定先路过它再遇到数据块；区域查询路径则在引擎初始化时按文件元信息预取。
     * 两条路径都保证工作线程取用时缓存已就位。
     */
    virtual AuxPayloadPtr getQualPrior(int64_t /*packageIndex*/) { return AuxPayloadPtr(); }

    /* 先验块容器头的绝对偏移；压缩侧据此写进块 meta，解压侧据此回查。-1 表示无先验。 */
    virtual int64_t getQualPriorAddress() const { return -1; }


    std::vector<AuxBlockConsumer*> auxConsumers;

    std::vector<std::thread> workThreads;
    std::thread writeThread;
    int64_t blockId2Write;

    uint32_t blockCount;

    /*
     * 块级失败的最终归口：工作线程里 actuatorProc 返回非零时置位，引擎收尾
     * （startEnginePostProc 之后）统一检查并返回错误，进程以非零退出——
     * 块失败不再是"打个警告、推个零长块就过去"。atomic 因为工作线程并发置位。
     */
    std::atomic<bool> taskFailed{false};

    /* 同步辅助块发射的交接槽位，由发射者与写线程共用，受 auxEmitMutex 保护。 */
    mutable std::mutex auxEmitMutex;
    mutable std::condition_variable auxEmitCond;
    int64_t auxEmitOffset = -1;
    bool auxEmitDone = false;

    /*
     * 首块串行化：id==0 的 coder 线程先处理完第 0 块再放行其余线程，避免并发
     * preAnalysis 在 SamInfo 尚未填充时把块误判为二进制。perf 分支原有机制，
     * 接入先验门闩（workStartBarrier）时被误删，此处恢复。放行只发生一次。
     */
    mutable std::mutex coderStartMutex;
    std::condition_variable coderStartCond;
    bool coderStartSync = false;

    /* 首次返回 true 并记下，之后恒返回 false；配合 notify_all 只放行一次。 */
    bool firstCoderNotify(bool flag);
};
