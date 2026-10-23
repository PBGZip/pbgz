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
#include "bam_columns.h"
#include "io_wrapper.h"
#include "pbgz_file_wrapper.h"
#include "utils/memory_util.h"
#include "bgzf_parallel.h"

namespace BlockUtil {
    bool isFastqBlock(BlockType type);
    bool isSAMBlock(BlockType type);
    bool isBAMBlock(BlockType type);

    /*
     * Whether a file is BAM: raw BAM magic, or a BGZF/gzip stream inflating to BAM.
     */
    bool isBamFile(const std::string& fileName);

    /*
     * Detect the input format from file content (SAM/BAM/FASTQ/BINARY/PBGZFILE), inflating gz/BGZF first. Regular files only; not for pipes/STDIN.
     */
    BlockType detectInputFileType(const std::string& fileName);

    /*
     * Auxiliary block: no user data, loaded on demand by seek; reference genome, its index and QUAL priors.
     */
    bool isAuxiliaryBlock(BlockType type);

    std::string getBlockTypeName(BlockType type);
}

/*
 * Block read base class: reads raw data from ioReader, then hands it to the subclass analyzeBlock for pre-analysis.
 */
class BlockReader {
public:
    /*
     * preReadData/preReadLen is the format-detection data prefetched by the Creator, merged into the first block.
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

    /* Read one data block: actual length after alignment, 0 for EOF, -1 for error. */
    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW);

    virtual int32_t init() { return 0; }

    /*
     * Per-block read target in bytes set by -l; 0 falls back to the block's own blockSize.
     */
    void setReadBlockBytes(uint32_t bytes) {
        readBlockBytes = bytes;
    }

    /*
     * Whether this block carries data lines; the SAM header block returns false.
     */
    virtual bool blockHasData(const RoughIOBlock* /*blockPtr*/) const { return true; }

protected:
    /* Per-block read target in bytes: prefer the -l target, otherwise the block's own blockSize */
    size_t readTargetBytes(const RoughIOBlock* blockPtr) const {
        return (readBlockBytes > 0) ? (size_t)readBlockBytes : blockPtr->getBlockSize();
    }

    /*
     * Pre-analyze the block for its format and return the block type; on failure Fastq/Sam return BINARY.
     */
    virtual BlockType analyzeBlock(RoughIOBlock* blockPtr, BlockType fileType);

    /*
     * Move the incomplete tail record into cache so coding lands on record boundaries; false on malformed tail.
     */
    bool alignToRecordBoundary(RoughIOBlock* blockPtr, bool isFastq);

    /* Merge prefetched data and leftover cache into the block buffer; returns bytes merged. */
    size_t prependBufferedData(RoughIOBlock* blockPtr, size_t& totalLen, bool capAtBlockSize = true);

protected:
    int64_t blockId;
    IOReader* ioReader;
    uint8_t cache[BLOCK_SIZE];  // Used to store remaining data from previous read
    uint64_t cacheLen;          // Length of data in cache
    uint8_t detectBuf[BLOCK_TYPE_DETECT_SIZE];  // Format-detection data prefetched by the Creator
    uint64_t detectLen;         // Length of the prefetched data
    uint32_t readBlockBytes = 0;   /* Per-block read target bytes from -l; 0 means use the block's own blockSize */
};

/* Binary (unrecognized) format */
class BinaryBlockReader : public BlockReader {
public:
    BinaryBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0)
        : BlockReader(reader, preReadData, preReadLen) { }

protected:
    virtual BlockType analyzeBlock(RoughIOBlock* /*blockPtr*/, BlockType /*fileType*/) override {
        return BINARY;
    }
};

/* FASTQ format */
class FastqBlockReader : public BlockReader {
public:
    FastqBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0)
        : BlockReader(reader, preReadData, preReadLen) { }

protected:
    virtual BlockType analyzeBlock(RoughIOBlock* blockPtr, BlockType fileType) override;
};


/* SAM format */
class SamBlockReader : public BlockReader {
public:
    /*
     * readsPerBlock: read-count bound per data block (1-5 -> 10000, 6-7 -> 25000, 8-9 -> 100000).
     * splitHeader=true gives the header its own block; false merges it with the first data block.
     */
    SamBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0,
                   uint32_t readsPerBlock = 10000, bool splitHeader = true)
        : BlockReader(reader, preReadData, preReadLen), readsPerBlock(readsPerBlock), splitHeader(splitHeader) {
        lastBlockHasData = false;
        /* CRAM rule: bases_per_slice = seqs_per_slice * 500. */
        blockMaxBases = (size_t)readsPerBlock * 500;
    }

    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW) override;

    virtual bool blockHasData(const RoughIOBlock* /*blockPtr*/) const override {
        return lastBlockHasData;
    }

protected:
    virtual BlockType analyzeBlock(RoughIOBlock* blockPtr, BlockType fileType) override;

protected:
    uint32_t readsPerBlock;      /* Maximum number of reads per data block */
    size_t blockMaxBases;        /* Base-count cap per data block (= readsPerBlock*500, CRAM rule) */
    bool splitHeader;            /* Whether the header forms its own block */
    bool lastBlockHasData;       /* Whether the most recently read block contains data lines (queried by the engine) */
};

/*
 * BAM format: the header becomes its own block (SAM text) and each alignment read is decompressed into a SAM line, grouped into readsPerBlock blocks. The block type stays BAM.
 */
class BamBlockReader : public BlockReader {
public:
    BamBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0,
                   uint32_t readsPerBlock = 10000, bool splitHeader = true)
        : BlockReader(reader, preReadData, preReadLen), readsPerBlock(readsPerBlock), splitHeader(splitHeader) {
        headerParsed = false;
        headerWritten = false;
        lastBlockHasData = false;
        /* CRAM rule: bases_per_slice = seqs_per_slice * 500. */
        blockMaxBases = (size_t)readsPerBlock * 500;
    }

    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW) override;

    virtual bool blockHasData(const RoughIOBlock* /*blockPtr*/) const override {
        return lastBlockHasData;
    }

    void setStructMode(bool on) { structMode = on; }

    /*
     * With rawStreamMode on, the reader renders the first textBlockLimit blocks as SAM text, then streams later blocks as raw records for the worker to render.
     */
    void setRawStreamMode(bool on, int64_t textBlockLimit) {
        rawStreamMode = on;
        rawStreamTextLimit = textBlockLimit;
    }

protected:
    /* Read the raw BAM byte stream (prefetched data + ioReader); BamGzBlockReader overrides this to inflate */
    virtual size_t readBamBytes(void* dst, size_t n);

    /* Read directly from prefetched data + ioReader (without inflate) */
    size_t readRawFromSource(void* dst, size_t n);

    int32_t parseBamHeader();
    int32_t parseBamRecord(const uint8_t* data, int32_t size, std::string& samLine);

    /*
     * Structured mode: append raw BAM record bytes only; BamCodecActuator parses them into columns on the worker thread.
     */
protected:
    std::vector<BamRef> refs;          /* BAM reference sequence list (name + length) */
    std::string headerText;            /* Generated SAM header text (including newlines) */
    bool headerParsed;
    bool headerWritten;
    uint32_t readsPerBlock;            /* Maximum number of reads per data block */
    size_t blockMaxBases;              /* Base-count cap per data block (= readsPerBlock*500, CRAM rule) */
    bool splitHeader;                  /* Whether the header forms its own block */
    bool lastBlockHasData;

    /*
     * A read whose SAM line would push the block past the budget is parked here for the next block.
     */
    std::string pendingSamLine;

    /* Reused scratch buffers for one BAM record and its parsed SAM line. */
    std::vector<uint8_t> recBuf;
    std::string lineBuf;
    bool structMode = false;
    std::shared_ptr<BamColumns> colsCur;      /* columns of the block being read */
    std::vector<uint8_t> pendingBamRecord;    /* a read parked to keep blocks inside the base cap */
    bool rawStreamMode = false;               /* archive: stream raw records after the first blocks (see setRawStreamMode) */
    int64_t rawStreamTextLimit = INT64_MAX;   /* blockId bound: blocks with id <= limit keep SAM-text rendering */
    std::shared_ptr<const std::vector<BamRef>> rawRefs; /* shared ref table handed to raw blocks */
};

/* GZ-compressed BAM: inflate the BGZF stream to raw BAM, then convert to SAM. */
class BamGzBlockReader : public BamBlockReader {
public:
    BamGzBlockReader(IOReader* reader, const uint8_t* preReadData = nullptr, uint64_t preReadLen = 0,
                     uint32_t readsPerBlock = 10000, bool splitHeader = true)
        : BamBlockReader(reader, preReadData, preReadLen, readsPerBlock, splitHeader) {
        inflateReady = false;
        gzInLen = 0;
        gzInEof = false;
        parStarted = false;
        parMode = false;
        parThreads = 4;    /* default until the caller supplies a budget (see setParThreads) */
    }

    ~BamGzBlockReader() {
        if (parStarted) {
            bgzfPar.stop();
        }
        if (inflateReady) {
            inflateEnd(&inflateState);
        }
    }

    /*
     * Number of inflate workers for the parallel BGZF decoder, derived from -t; a missing value keeps the default of 4.
     */
    void setParThreads(int n) { parThreads = (n >= 1) ? n : 4; }

protected:
    virtual size_t readBamBytes(void* dst, size_t n) override;

    /* Serial fallback used when the stream is not cuttable BGZF. */
    size_t readRawBam(void* dst, size_t n);

private:
    z_stream inflateState;
    bool inflateReady;
    uint8_t gzInBuf[1 << 20];          /* BGZF input buffer (1MB: fewer read/refill rounds) */
    size_t gzInLen;
    bool gzInEof;

    BgzfParallelDecoder bgzfPar;
    bool parStarted;                   /* probe done */
    bool parMode;                      /* running the parallel decoder */
    int parThreads;
    std::vector<uint8_t> probeBuf;     /* stream head read once for the probe */
};

class BlockWriter;
class PbgzBlockWriter;

/*
 * BGZF inflate workers on the compression read path: one per six coding threads,
 * rounded up (at least one).
 */
inline uint32_t bamInflateThreadNum(uint32_t threadNum)
{
    return (threadNum > 0) ? ((threadNum + 5) / 6) : 1;
}

/*
 * BGZF deflate workers on the decompression write path: one per three coding
 * threads, rounded up (at least one).
 */
inline uint32_t bamDeflateThreadNum(uint32_t threadNum)
{
    return (threadNum > 0) ? ((threadNum + 2) / 3) : 1;
}

/*
 * BAM<->SAM conversion helper budget (BAM-read render stage and the -b sam2bam
 * stage): one per eight coding threads, rounded up (at least one).
 */
inline uint32_t bamConvertThreadNum(uint32_t threadNum)
{
    return (threadNum > 0) ? ((threadNum + 7) / 8) : 1;
}

/*
 * Creator for BlockReader/BlockWriter: reads a small sample to pick the format, then creates the matching reader; the block type is always the true format, never a GZ variant.
 */
class BlockFactory {
public:
    /*
     * compressLevel sets SAM block granularity (1-5 -> 10000, 6-7 -> 25000, 8-9 -> 100000 reads).
     * splitSamHeader=false merges the header into the first data block. bamStructMode reads raw records for column parsing. pbgzAsOpaque reads an existing pbgz archive as opaque binary. bgzfThreads sets the BGZF inflate worker count (0 keeps the default).
     */
    static BlockReader* createBlockReader(IOReader* ioReader, uint8_t compressLevel = 0,
                                          bool splitSamHeader = true,
                                          bool bamStructMode = false,
                                          bool pbgzAsOpaque = false,
                                          uint32_t bgzfThreads = 0);

    /*
     * The SAM block granularity used for a level (1-5 -> 10000, 6-7 -> 25000, 8-9 -> 100000), so pre-selection trials at real block volume.
     */
    static uint32_t samReadsPerBlockOfLevel(uint8_t compressLevel);

    static BlockWriter* createBlockWriter(IOWriter* ioWriter);

    static PbgzBlockWriter* createPbgzBlockWriter(IOWriter* ioWriter);
};

class PbgzBlockReader : public BlockReader {
public:
    /*
     * preReadData/preReadLen is the prefetched format-detection data, consumed by the internal PbgzFileReader.
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

    // Sticky error from the underlying writes, retrieved by the writer thread before exiting and handed to the engine for final handling
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
 * Convert a decompressed SAM text block into BAM records in place, off the writer thread.
 * The first block also emits the BAM header; blocks already carrying BAM bytes (structured
 * BAM) are left untouched. Single-threaded by contract: one converter instance is used by
 * exactly one stage thread.
 */
class SamToBamConverter {
public:
    explicit SamToBamConverter(uint32_t maxThreads = 1) : maxThreads(maxThreads) { }

    /* Returns 0 on success. On success the block is marked bamPrebuilt when it now holds BAM bytes. */
    int32_t convert(RoughIOBlock* blockPtr);

private:
    bool headerWritten = false;
    bool passThrough = false;              /* First block was not SAM: leave the stream untouched. */
    std::map<std::string, int32_t> refIndex;
    uint32_t maxThreads;                   /* Data-line conversion budget inside a block. */
};

/*
 * BAM output writer (-b): frames BAM byte blocks into BGZF and writes them.
 * With the -b conversion stage blocks arrive prebuilt; text blocks are still handled
 * here as a fallback (the first block yields the BAM header).
 */
class BamWriter : public BlockWriter {
public:
    /*
     * maxThreads caps the writer's convert and deflate threads; 0 means hardware_concurrency.
     */
    explicit BamWriter(IOWriter* pIoWriter, uint32_t maxThreads = 0);

    virtual ~BamWriter();

    virtual int32_t writeBlock(RoughIOBlock* blockPtr) override;

    /*
     * Flush the BGZF buffer and append the EOF marker; idempotent.
     */
    int32_t finish();

private:
    /* Parse the reference sequence list from SAM header lines and write the BAM header; returns pos (offset of the first data line) or dataLen */
    int32_t writeBamHeader(const uint8_t* data, size_t len, size_t& dataStart);

    /* Process the data lines in a span of SAM text (converting each into a BAM record) */
    int32_t writeDataLines(const uint8_t* buffer, size_t start, size_t end);

    /* Write as-is (pass-through mode, or the underlying write for BGZF) */
    int32_t writeRaw(const void* data, size_t len);

    /* Buffer uncompressed BAM bytes and frame every complete BGZF block */
    int32_t bgzfWrite(const void* data, size_t len);

    /*
     * Frame the buffered bytes into BGZF blocks, deflating complete 64KB blocks in parallel.
     * bgzfWrite() emits every complete block; bgzfFlushRemaining() emits the trailing partial one.
     */
    int32_t bgzfFlushRemaining();

    /* Sink called by the deflate pool (on the writer thread) for each completed BGZF block. */
    static int32_t bgzfSink(void* ctx, const uint8_t* data, size_t len);

private:
    std::map<std::string, int32_t> refIndex;  /* Reference name -> refID */

    bool headerWritten;
    bool passThrough;                      /* Set when the first block is not SAM; pass through as-is */
    bool finished;                         /* finish() has been executed (prevents duplicate EOF markers) */

    /*
     * Uncompressed BAM bytes not yet framed into BGZF blocks; kept across writeBlock
     * calls when the deflate pool could not be started (serial fallback).
     */
    std::vector<uint8_t> bgzfPending;

    /* Independent BGZF deflate worker pool; the writer only submits/copies and writes results. */
    BgzfParallelEncoder bgzfEncoder;
    bool bgzfEncoderReady;

    uint32_t maxThreads;   /* 0 = hardware_concurrency (see the constructor) */
};

/*
 * Convert one BAM alignment record into a SAM line; `refs` resolves RNAME/RNEXT. Returns non-zero on error.
 */
int32_t bamRecordToSamLine(const uint8_t* data, int32_t size,
                           const std::vector<BamRef>& refs, std::string& line);

/*
 * Worker-side render of a raw-record BAM block: rebuild the SAM text, recompute npos and set the MD5. Returns 0 on success.
 */
int32_t bamRenderRawBlockToText(RoughIOBlock* blockPtr);
