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
     * Determine whether a file is BAM: a raw BAM magic number, or a BGZF/gzip stream that
     * inflates to BAM.
     *
     * A standard BAM is a BGZF (gzip) stream and is byte-for-byte indistinguishable from a .gz
     * file; the engine uses this to decide whether input goes through "transparent gz
     * decompression" or is handed to FileReader + BamGzBlockReader. BAM's gzip is inside the
     * format; if transparent decompression were used, FileReader would not know the total file
     * size, and the "is it worth writing" heuristic for the QUAL prior would degrade.
     */
    bool isBamFile(const std::string& fileName);

    /*
     * Detect the input format from file content (SAM/BAM/FASTQ/BINARY/PBGZFILE), auto-inflating
     * gz/BGZF before judging. Used for pre-compression decisions (e.g. loading a reference for
     * SAM/BAM only needs squash, no read-mapping hash table). Detection only reads the file
     * header and does not modify the file; callers must not use it on pipes/STDIN.
     */
    BlockType detectInputFileType(const std::string& fileName);

    /*
     * Auxiliary block: carries no user data and does not enter the "read -> parallel compress ->
     * write in order" pipeline; it is loaded on demand by seek using file-level offsets (reference
     * genome, its index, and QUAL priors all belong here). Pipeline entry points only query this
     * single predicate, so adding a new auxiliary block type requires no new branches at each entry.
     */
    bool isAuxiliaryBlock(BlockType type);

    std::string getBlockTypeName(BlockType type);
}

/*
 * Block read base class: reads raw data from ioReader (including the format-detection data
 * prefetched by the Creator and the cache left by the previous block read), then hands it to
 * the subclass analyzeBlock for format pre-analysis.
 *
 * There is one subclass per file format (Fastq/Sam/Bam/FastqGz/BamGz/Binary/Pbgz); the concrete
 * type is determined and created by BlockFactory after prefetching.
 */
class BlockReader {
public:
    /*
     * preReadData/preReadLen is the format-detection data the Creator read before creating this
     * reader; on the first readBlock call it is merged verbatim into the target RoughIOBlock's buffer.
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
     * Read one data block. Returns the actual data length (after alignment); 0 means EOF, -1 means a read error.
     */
    virtual int64_t readBlock(RoughIOBlock* blockPtr, BlockType fileType = TYPE_UNKNOW);

    virtual int32_t init() { return 0; }

    /*
     * Set the per-block read target in bytes determined by the -l option (the block initially
     * allocates only 1MB and grows to this target on demand). 0 means falling back to the
     * block's own blockSize.
     */
    void setReadBlockBytes(uint32_t bytes) {
        readBlockBytes = bytes;
    }

    /*
     * Whether this block carries data lines. Default true; the SAM header block (only @ lines)
     * returns false, which lets the engine defer file-level decisions (codec selection) until
     * the first data block.
     */
    virtual bool blockHasData(const RoughIOBlock* /*blockPtr*/) const { return true; }

protected:
    /* Per-block read target in bytes: prefer the -l target, otherwise the block's own blockSize */
    size_t readTargetBytes(const RoughIOBlock* blockPtr) const {
        return (readBlockBytes > 0) ? (size_t)readBlockBytes : blockPtr->getBlockSize();
    }

    /*
     * Subclasses pre-analyze the block content according to their format (newline position
     * recording, Fastq/Sam actuator preAnalysis) and return the determined block type. On
     * pre-analysis failure, Fastq/Sam must return BINARY so that compression creates a binary
     * Actuator.
     */
    virtual BlockType analyzeBlock(RoughIOBlock* blockPtr, BlockType fileType);

    /*
     * Move the incomplete record at the block tail into cache so that structured coding lands
     * exactly on record boundaries. Returns false when the tail block (dataLen < read target)
     * contains malformed data; the caller should degrade to BINARY.
     */
    bool alignToRecordBoundary(RoughIOBlock* blockPtr, bool isFastq);

    /* Merge the Creator's prefetched data and the leftover cache into the block buffer and
     * return the number of bytes merged. With capAtBlockSize=false (SAM blocks split by read
     * line count may exceed byte blockSize), all prefetched data and cache are merged at once,
     * growing the buffer on demand. */
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
     * readsPerBlock: upper bound for splitting the data region by read line count (set by the
     * compression level: 1-5 -> 10000, 6-7 -> 25000, 8-9 -> 100000).
     * splitHeader=true: header lines form their own block (blockId==0 returns only @ lines);
     * =false (downstream like sorting needs self-contained @SQ within a block) merges the header
     * with the first data block, still splitting the data region by readsPerBlock.
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
    uint32_t readsPerBlock;      /* Maximum number of reads per data block */
    bool splitHeader;            /* Whether the header forms its own block */
    bool lastBlockHasData;       /* Whether the most recently read block contains data lines (queried by the engine) */
};

/*
 * BAM format: reads the raw BAM byte stream, returns the header as an independent block
 * (converted to SAM header text), decompresses each alignment read into a SAM line, groups
 * readsPerBlock reads into a SAM block, and hands it to the SAM compressor. The block type is
 * marked BAM (FileType prints BAM), but the content is SAM text, so every capability SAM
 * supports (codec, priors, indexing, etc.) is also supported for BAM.
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
    /* Read the raw BAM byte stream (prefetched data + ioReader); BamGzBlockReader overrides this to inflate */
    virtual size_t readBamBytes(void* dst, size_t n);

    /* Read directly from prefetched data + ioReader (without inflate) */
    size_t readRawFromSource(void* dst, size_t n);

    int32_t parseBamHeader();
    int32_t parseBamRecord(const uint8_t* data, int32_t size, std::string& samLine);

protected:
    struct BamRef { std::string name; int32_t len; };
    std::vector<BamRef> refs;          /* BAM reference sequence list (name + length) */
    std::string headerText;            /* Generated SAM header text (including newlines) */
    bool headerParsed;
    bool headerWritten;
    uint32_t readsPerBlock;            /* Maximum number of reads per data block */
    bool splitHeader;                  /* Whether the header forms its own block */
    bool lastBlockHasData;
};

/* GZ-compressed BAM format: the inner layer is still a BGZF stream; inflate it to raw BAM first, then convert to SAM */
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
    uint8_t gzInBuf[64 * 1024];        /* BGZF input buffer */
    size_t gzInLen;
    bool gzInEof;
};

class BlockWriter;
class PbgzBlockWriter;

/*
 * Creator for BlockReader/BlockWriter.
 *
 * createBlockReader first reads a small amount of data to determine the file format
 * (FASTQ/SAM/BAM/PBGZ, treating unrecognizable input as binary), then creates the matching
 * BlockReader subclass; prefetched data is handed to the subclass to merge into the first block.
 * GZ-compressed input is transparently decompressed by the io layer or inflated inside the
 * reader; the block type is always set to the true format after decompression
 * (FASTQ/SAM/BAM/BINARY), never a GZ variant.
 */
class BlockFactory {
public:
    /*
     * compressLevel determines the read-line-count granularity for splitting SAM/SAM-GZ data
     * blocks: 1-5 -> 10000 reads/block, 6-7 -> 25000, 8-9 -> 100000. With splitSamHeader=false
     * the SAM header does not form its own block (merging with the first data block), for
     * downstream consumers such as sorting that need a self-contained @SQ within the block.
     */
    static BlockReader* createBlockReader(IOReader* ioReader, uint8_t compressLevel = 0,
                                          bool splitSamHeader = true);

    static BlockWriter* createBlockWriter(IOWriter* ioWriter);

    static PbgzBlockWriter* createPbgzBlockWriter(IOWriter* ioWriter);
};

class PbgzBlockReader : public BlockReader {
public:
    /*
     * preReadData/preReadLen is the format-detection data prefetched by the Creator (piped
     * input cannot be seeked back, so it must be handed over); the internal PbgzFileReader
     * consumes these bytes before reading from ioReader.
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
 * BAM output writer: converts decompressed SAM text blocks into standard BAM (BGZF-compressed)
 * and writes them to the underlying ioWriter.
 *
 * Usage: DecompressEngine selects -b to use it in place of the plain BlockWriter when
 * decompressing.
 *   - The first SAM block (@ header lines, possibly with merged data lines) is parsed to extract
 *     the reference sequence list and writes the BAM header;
 *   - Subsequent data blocks are converted line by line into BAM records;
 *   - Output is compressed in BGZF blocks; at file end the residual block is flushed and a BGZF
 *     EOF marker is appended (standard BAM).
 *
 * If the first block is not SAM/BAM (e.g. decompressing FASTQ/binary), it is not converted but
 * passed through as-is with a warning, avoiding data loss.
 */
class BamWriter : public BlockWriter {
public:
    explicit BamWriter(IOWriter* pIoWriter);

    virtual ~BamWriter();

    virtual int32_t writeBlock(RoughIOBlock* blockPtr) override;

    /*
     * Finish: compress the remaining data in the BGZF buffer into a block and append the BGZF
     * EOF marker. Called before the writer thread exits (releaseBlockWriter); the destructor
     * calls it once more as a fallback (idempotent).
     */
    int32_t finish();

private:
    /* Parse the reference sequence list from SAM header lines and write the BAM header; returns pos (offset of the first data line) or dataLen */
    int32_t writeBamHeader(const uint8_t* data, size_t len, size_t& dataStart);

    /* Convert one SAM alignment line into a BAM record and write it */
    int32_t writeBamRecord(const uint8_t* line, size_t len);

    /* Process the data lines in a span of SAM text (converting each into a BAM record) */
    int32_t writeDataLines(const uint8_t* buffer, size_t start, size_t end);

    /* Write as-is (pass-through mode, or the underlying write for BGZF) */
    int32_t writeRaw(const void* data, size_t len);

    /* Write through BGZF block compression (the container of standard BAM) */
    int32_t bgzfWrite(const void* data, size_t len);

    /* Compress the bytes accumulated in the buffer into one BGZF block and write it */
    int32_t bgzfFlushBlock();

private:
    struct BamRef {
        std::string name;
        int32_t len;
    };

    std::vector<BamRef> refs;              /* BAM reference sequence list (matching @SQ) */
    std::map<std::string, int32_t> refIndex;  /* Reference name -> refID */

    bool headerWritten;
    bool passThrough;                      /* Set when the first block is not SAM; pass through as-is */
    bool finished;                         /* finish() has been executed (prevents duplicate EOF markers) */

    /* BGZF block compression state */
    uint8_t bgzfBuf[65536];
    size_t bgzfLen;
    z_stream bgzfZs;
    bool bgzfReady;
};
