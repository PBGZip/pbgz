/*
 * block_wrapper.h - Header file for pbgz project
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

#include <cstring>
#include <map>

#include "io_block.h"
#include "io_wrapper.h"
#include "pbgz_file_wrapper.h"
#include "utils/memory_util.h"

namespace BlockUtil {
    bool isFastqBlock(BlockType type);
    bool isSAMBlock(BlockType type);
    bool isBAMBlock(BlockType type);

    /*
     * 判定一个文件是否是 BAM：原始 BAM 魔数，或 BGZF/gzip 压缩后解压即 BAM。
     *
     * 标准 BAM 是 BGZF（gzip）流，逐字节看和 .gz 文件无法区分；引擎据此决定输入是
     * 走"透明 gz 解压"还是交给 FileReader + BamGzBlockReader。BAM 的 gzip 是格式内层，
     * 若走了透明解压，FileReader 拿不到文件总长，QUAL 先验的"值不值得写"判据会退化。
     */
    bool isBamFile(const std::string& fileName);

    /*
     * 按文件内容探测输入格式（SAM/BAM/FASTQ/BINARY/PBGZFILE），gz/BGZF 自动解压探测。
     * 用于压缩前决策（如 SAM/BAM 加载参考只需 squash，无需建读映射哈希表）。
     * 探测只读文件头部，不改动文件；调用方不得对管道/STDIN 使用。
     */
    BlockType detectInputFileType(const std::string& fileName);

    /*
     * 辅助块：不携带用户数据，不进"读→并行压缩→按序写"这条流水线，
     * 由文件级偏移按需 seek 加载（参考基因组、其索引、QUAL 先验都属此类）。
     * 流水线入口只问这一个谓词，新增辅助块类型不必再去每个入口补分支。
     */
    bool isAuxiliaryBlock(BlockType type);

    std::string getBlockTypeName(BlockType type);
}

/*
 * 块读取基类：负责从 ioReader 读原始数据（含 Creator 预读的格式探测数据、
 * 上次读块留下的缓存），再交给子类 analyzeBlock 完成格式预分析。
 *
 * 每个文件格式一个子类（Fastq/Sam/Bam/FastqGz/BamGz/Binary/Pbgz），
 * 具体类型由 BlockFactory 在预读后判定并创建。
 */
class BlockReader {
public:
    /*
     * preReadData/preReadLen 是 Creator 创建本 reader 前读出的格式探测数据，
     * readBlock 首次调用时会原样并入目标 RoughIOBlock 的缓冲。
     */
    BlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0) {
        ioReader = reader;
        blockId = 0;
        cacheLen = 0;
        memset(cache, 0, sizeof(cache));
        detectLen = 0;
        if (preReadData != nullptr && preReadLen > 0 && preReadLen <= BLOCK_TYPE_DETECT_SIZE) {
            memcpy(detectBuf, preReadData, preReadLen);
            detectLen = preReadLen;
        }
    }

    virtual ~BlockReader() {
        ioReader = nullptr;
    }

    /*
     * 读一个数据块。返回实际数据长度（对齐后），0 表示 EOF，-1 表示读错误。
     */
    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW);

    virtual int32_t init() { return 0; }

    /*
     * 设置按 -l 参数决定的单块读取目标字节数（块初始只分配 1MB，读取时按此目标
     * 按需扩容）。0 表示回落到块自身的 blockSize。
     */
    void setReadBlockBytes(uint32_t bytes) {
        readBlockBytes = bytes;
    }

    /*
     * 本块是否携带数据行。默认 true；SAM 头部块（只含 @ 行）返回 false，
     * 引擎据此把文件级决策（编码器选型）推迟到第一个数据块上执行。
     */
    virtual bool blockHasData(const RoughIOBlock* /*blockPtr*/) const { return true; }

protected:
    /* 单块读取目标字节数：优先用 -l 决定的目标，否则用块自身 blockSize */
    size_t readTargetBytes(const RoughIOBlock* blockPtr) const {
        return (readBlockBytes > 0) ? (size_t)readBlockBytes : blockPtr->getBlockSize();
    }

    /*
     * 子类按各自格式对块内容做预分析（换行符位置记录、Fastq/Sam actuator 的
     * preAnalysis），返回判定出的块类型。Fastq/Sam 预分析失败时须返回 BINARY，
     * 便于后续压缩创建二进制 Actuator。
     */
    virtual BlockType analyzeBlock(RoughIOBlock* blockPtr, BlockType fileType);

    /*
     * 把块尾不完整的记录移入 cache，保证结构化编码恰好落在记录边界上。
     * 末尾块（dataLen < 读取目标）数据不合规时返回 false，调用方应降级为 BINARY。
     */
    bool alignToRecordBoundary(RoughIOBlock* blockPtr, bool isFastq);

    /* 把 Creator 预读数据与上次遗留缓存并入块缓冲，返回并入字节数。
     * capAtBlockSize=false 时（SAM 按 read 行数分块，块可超过 byte blockSize）
     * 会一次性并入全部预读数据与缓存，并按需扩容。 */
    size_t prependBufferedData(RoughIOBlock* blockPtr, size_t& totalLen, bool capAtBlockSize = true);

protected:
    int64_t blockId;
    IOReader* ioReader;
    uint8_t cache[BLOCK_SIZE];  // Used to store remaining data from previous read
    uint64_t cacheLen;          // Length of data in cache
    uint8_t detectBuf[BLOCK_TYPE_DETECT_SIZE];  // Creator 预读的格式探测数据
    uint64_t detectLen;         // 预读数据长度
    uint32_t readBlockBytes = 0;   /* -l 决定的单块读取目标字节数；0 表示用块自身 blockSize */
};

/* 二进制（无法识别）格式 */
class BinaryBlockReader : public BlockReader {
public:
    BinaryBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0)
        : BlockReader(reader, preReadData, preReadLen) { }

protected:
    virtual BlockType analyzeBlock(RoughIOBlock* /*blockPtr*/, BlockType /*fileType*/) override {
        return BINARY;
    }
};

/* FASTQ 格式 */
class FastqBlockReader : public BlockReader {
public:
    FastqBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0)
        : BlockReader(reader, preReadData, preReadLen) { }

protected:
    virtual BlockType analyzeBlock(RoughIOBlock* blockPtr, BlockType fileType) override;
};


/* SAM 格式 */
class SamBlockReader : public BlockReader {
public:
    /*
     * readsPerBlock：数据区按 read 行数分块的上界（由压缩级别决定：1-5 -> 10000，
     * 6-7 -> 25000，8-9 -> 100000）。
     * splitHeader=true：头部行独立成块（blockId==0 只返回 @ 行）；=false（排序等下游
     * 需要块内自含 @SQ）时头部与首个数据块合并，但仍按 readsPerBlock 切分数据区。
     */
    SamBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0,
                   uint32_t readsPerBlock = 10000, bool splitHeader = true)
        : BlockReader(reader, preReadData, preReadLen), readsPerBlock(readsPerBlock), splitHeader(splitHeader) {
        lastBlockHasData = false;
    }

    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW) override;

    virtual bool blockHasData(const RoughIOBlock* /*blockPtr*/) const override {
        return lastBlockHasData;
    }

protected:
    virtual BlockType analyzeBlock(RoughIOBlock* blockPtr, BlockType fileType) override;

protected:
    uint32_t readsPerBlock;      /* 每个数据块最多包含的 read 数 */
    bool splitHeader;            /* 头部是否独立成块 */
    bool lastBlockHasData;       /* 最近一次读出的块是否含数据行（引擎问询用） */
};

/*
 * BAM 格式：读出原始 BAM 字节流，头部作为独立块返回（转成 SAM 头文本），
 * 比对区逐条 read 解压成 SAM 行，按 readsPerBlock 条 read 组成一个 SAM 块，
 * 交给 SAM 压缩器压缩。块类型标记为 BAM（打印 FileType 显示 BAM），但内容
 * 是 SAM 文本，因此一切 SAM 支持的能力（编解码、先验、索引等）BAM 都支持。
 */
class BamBlockReader : public BlockReader {
public:
    BamBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0,
                   uint32_t readsPerBlock = 10000, bool splitHeader = true)
        : BlockReader(reader, preReadData, preReadLen), readsPerBlock(readsPerBlock), splitHeader(splitHeader) {
        headerParsed = false;
        headerWritten = false;
        lastBlockHasData = false;
    }

    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW) override;

    virtual bool blockHasData(const RoughIOBlock* /*blockPtr*/) const override {
        return lastBlockHasData;
    }

protected:
    /* 读取原始 BAM 字节流（预读数据 + ioReader）；BamGzBlockReader 覆写为 inflate */
    virtual size_t readBamBytes(void* dst, size_t n);

    /* 直接从预读数据 + ioReader 读（不经过 inflate） */
    size_t readRawFromSource(void* dst, size_t n);

    int32_t parseBamHeader();
    int32_t parseBamRecord(const uint8_t* data, int32_t size, std::string& samLine);

protected:
    struct BamRef { std::string name; int32_t len; };
    std::vector<BamRef> refs;          /* BAM 参考序列表（name + length） */
    std::string headerText;            /* 生成的 SAM 头部文本（含换行） */
    bool headerParsed;
    bool headerWritten;
    uint32_t readsPerBlock;            /* 每个数据块最多包含的 read 数 */
    bool splitHeader;                  /* 头部是否独立成块 */
    bool lastBlockHasData;
};

/* GZ 压缩的 BAM 格式：内层仍是 BGZF 流，先 inflate 成原始 BAM 再转 SAM */
class BamGzBlockReader : public BamBlockReader {
public:
    BamGzBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0,
                     uint32_t readsPerBlock = 10000, bool splitHeader = true)
        : BamBlockReader(reader, preReadData, preReadLen, readsPerBlock, splitHeader) {
        inflateReady = false;
        gzInLen = 0;
        gzInEof = false;
    }

    ~BamGzBlockReader() {
        if (inflateReady) {
            inflateEnd(&inflateState);
        }
    }

protected:
    virtual size_t readBamBytes(void* dst, size_t n) override;

private:
    z_stream inflateState;
    bool inflateReady;
    uint8_t gzInBuf[64 * 1024];        /* BGZF 输入缓冲 */
    size_t gzInLen;
    bool gzInEof;
};

class BlockWriter;
class PbgzBlockWriter;

/*
 * BlockReader/BlockWriter 创建器。
 *
 * createBlockReader 先读一小部分数据解析文件格式（FASTQ/SAM/BAM/PBGZ，
 * 无法识别则视为二进制），再按格式创建对应的 BlockReader 子类；预读数据会
 * 交给子类并入首个块。GZ 压缩的输入由 io 层透明解压或 reader 内部 inflate，
 * 块类型一律按解压后的真实格式设置（FASTQ/SAM/BAM/BINARY），不设 GZ 变体。
 */
class BlockFactory {
public:
    /*
     * compressLevel 决定 SAM/SAM-GZ 数据块按 read 行数的分块粒度：
     * 1-5 -> 10000 read/块，6-7 -> 25000，8-9 -> 100000。
     * splitSamHeader=false 时 SAM 头部不独立成块（与首个数据块合并），
     * 供排序等需要块内自含 @SQ 的下游使用。
     */
    static BlockReader* createBlockReader(IOReader* ioReader, uint8_t compressLevel = 0,
                                          bool splitSamHeader = true);

    static BlockWriter* createBlockWriter(IOWriter* ioWriter);

    static PbgzBlockWriter* createPbgzBlockWriter(IOWriter* ioWriter);
};

class PbgzBlockReader : public BlockReader {
public:
    /*
     * preReadData/preReadLen 是 Creator 预读的格式探测数据（管道输入不可 seek 回退，
     * 必须交还）；内部 PbgzFileReader 会先消费这些字节再读 ioReader。
     */
    PbgzBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0)
        : BlockReader(reader) {
        pbgzFileReader = MemoryUtil::safeNewClass<PbgzFileReader>(ioReader);
        pbgzFileReader->setPreReadData(preReadData, preReadLen);
    }

    ~PbgzBlockReader() {
       MemoryUtil::safeDeleteClass(pbgzFileReader);
    }

    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW) override;

    virtual int32_t init() override;

    PbgzFileMeta& getBaseFileMeta() {
        return pbgzFileReader->getBaseFileMeta();
    }

    PbgzFileMeta& getDynamicFileMeta() {
        return pbgzFileReader->getDynamicFileMeta();
    }

    PbgzFileHeader& getFileHeader() {
         return pbgzFileReader->getFileHeader();
    }

    uint64_t getCurrentFileStart() const {
        return pbgzFileReader->getCurrentFileStart();
    }

    uint64_t getCurrentBlockStart() const {
        return pbgzFileReader->getCurrentBlockStart();
    }

    int32_t getCurrentFileIndex() const {
        return pbgzFileReader->getCurrentFileIndex();
    }

private:
    PbgzFileReader* pbgzFileReader;
};


class BlockWriter {
public:
    BlockWriter(IOWriter* pIoWriter): ioWriter(pIoWriter) { }

    virtual int32_t writeBlock(RoughIOBlock* blockPtr);

    virtual int32_t init() { return 0; }

    virtual ~BlockWriter() {
        ioWriter = nullptr;
    }

    // 底层写出的粘性错误，供写线程退出前取走交给引擎做最终处理
    int32_t getWriteError() const { return ioWriter != nullptr ? ioWriter->getWriteError() : 0; }

protected:
    IOWriter* ioWriter;
};


class PbgzBlockWriter : public BlockWriter {
public:
    PbgzBlockWriter(IOWriter* pIoWriter) : BlockWriter(pIoWriter) {
        pbgzFileWriter = MemoryUtil::safeNewClass<PbgzFileWriter>(ioWriter);;
    }
    virtual int32_t writeBlock(RoughIOBlock* blockPtr);

    void setBaseFileMeta(PbgzFileMeta& fileMeta) {
        pbgzFileWriter->setBaseFileMeta(fileMeta);
    }

    void setDynamicFileMeta(PbgzFileMeta& fileMeta) {
        pbgzFileWriter->setDynamicFileMeta(fileMeta);
    }

    virtual int32_t init();

    virtual ~PbgzBlockWriter() {
        MemoryUtil::safeDeleteClass(pbgzFileWriter);
    }

    int32_t writeBaseFileMeta();

    int32_t writeDynamicFileMeta();

    void updateHeadExt();

private:
    PbgzFileWriter* pbgzFileWriter;
};

/*
 * BAM 输出写入器：把解压出的 SAM 文本块转换成标准 BAM（BGZF 压缩）写入底层 ioWriter。
 *
 * 用法：DecompressEngine 在解压时指定 -b 用它替代普通 BlockWriter。
 *   - 首个 SAM 块（@ 头部行，可含合并的数据行）解析出参考序列表并写出 BAM 头；
 *   - 后续数据块逐行转换成 BAM 记录；
 *   - 输出按 BGZF 块压缩，文件收尾时 flush 残留块并追加 BGZF EOF 标记（标准 BAM）。
 *
 * 若首个块不是 SAM/BAM（如解压的是 FASTQ/二进制），则不转换、原样透传并告警，
 * 避免数据丢失。
 */
class BamWriter : public BlockWriter {
public:
    explicit BamWriter(IOWriter* pIoWriter);

    virtual ~BamWriter();

    virtual int32_t writeBlock(RoughIOBlock* blockPtr) override;

    /*
     * 收尾：把 BGZF 缓冲里剩余的数据压成块写出，并追加 BGZF EOF 标记。
     * 写线程退出前（releaseBlockWriter）调用；析构时兜底再调一次（幂等）。
     */
    int32_t finish();

private:
    /* 从 SAM 头部行解析参考序列表并写出 BAM 头；返回 pos（首个数据行的偏移）或 dataLen */
    int32_t writeBamHeader(const uint8_t* data, size_t len, size_t& dataStart);

    /* 把一条 SAM 比对行转成 BAM 记录写出 */
    int32_t writeBamRecord(const uint8_t* line, size_t len);

    /* 处理一段 SAM 文本里的数据行（逐行转成 BAM 记录） */
    int32_t writeDataLines(const uint8_t* buffer, size_t start, size_t end);

    /* 原样写出（透传模式，或 BGZF 的底层写） */
    int32_t writeRaw(const void* data, size_t len);

    /* 经 BGZF 块压缩写出（标准 BAM 的容器） */
    int32_t bgzfWrite(const void* data, size_t len);

    /* 把缓冲里积攒的字节压成一个 BGZF 块写出 */
    int32_t bgzfFlushBlock();

private:
    struct BamRef {
        std::string name;
        int32_t len;
    };

    std::vector<BamRef> refs;              /* BAM 参考序列表（与 @SQ 一致） */
    std::map<std::string, int32_t> refIndex;  /* 参考名 -> refID */

    bool headerWritten;
    bool passThrough;                      /* 首个块不是 SAM 时置真，原样透传 */
    bool finished;                         /* finish() 已执行（防重复 EOF 标记） */

    /* BGZF 块压缩状态 */
    uint8_t bgzfBuf[65536];
    size_t bgzfLen;
    z_stream bgzfZs;
    bool bgzfReady;
};
