/*
 * reference.h - Header file for pbgz project
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

#include <stdint.h>
#include <string>
#include <map>
#include <memory>
#include <json/json.h>
#include "utils/guard_bar.h"
#include "io_wrapper.h"

typedef std::vector<std::pair<std::pair<uint32_t*, uint32_t>, uint32_t>> HashTable;
typedef std::pair<std::pair<int64_t, uint8_t*>, std::pair<int64_t, int64_t>> RefeInfo;

class Reference {
private:
    typedef struct {
        uint32_t baseGroupPos;
        uint32_t hashBucket;
    } BaseGroupHash;

public:
    Reference(const std::string& fasqaName, uint32_t threadNum);

    virtual ~Reference();

    /// @brief Create index table
    /// @return  bool true for success / false for failure
    bool makeIndex();

    /// @brief Query position information in reference genome based on hash value
    /// @param hash Input hash value
    /// @param length Output parameter, returns the number of position information in this hash bucket
    /// @return Pointer to position information array, array contains length uint32_t values
    const uint32_t* queryPosition(const uint32_t &hash, uint32_t &length);

    /// @brief Get base group length used when building index
    /// @return Base group length, fixed to 31
    uint32_t getBaseGroupLength() const;

    /// @brief Get base group step used when building index
    /// @return Base group step, fixed to 32
    uint32_t getBaseGroupStep() const;

    /// @brief Get compressed data buffer of reference genome
    /// @return Pointer to compressed data buffer, each byte contains 2-bit encoding of 4 bases
    const uint8_t* getSquash() const;

    /// @brief Get length of reference genome compressed data
    /// @return Byte length of compressed data
    int64_t getSquashLength() const;

    /// @brief Get FASTA reference genome file name
    /// @return String reference of FASTA file name
    const std::string& getFastaFileName() const;

    /// @brief Get content length of FASTA reference genome file
    /// @return Byte length of FASTA file
    int64_t getFastaLength() const;

    /// @brief Get MD5 checksum of FASTA reference genome file
    /// @return MD5 checksum string
    const std::string& getFastaChecksum() const;

    /// @brief Clean up unmatched regions in reference genome
    /// Set unmatched squash bytes to 0 within specified range
    /// @param startPos Starting squash position (byte index)
    /// @param length Length to be cleaned (in bytes)
    void sanitizeRefSquash(int64_t startPos, int64_t length);

    /// @brief Update matched gene region information
    /// Mark specified position range of reference genome as matched
    /// @param actgPos Starting ACTG base position
    /// @param matchLength Matched base length
    void updateMatchedGene(uint64_t actgPos, uint32_t matchLength);

    /// @brief Build a NI index file from the FASTA reference at an explicit path
    /// Only ever invoked on explicit user request; nothing is cached implicitly.
    /// @param niFile NI file path to create
    /// @return true for success, false for failure
    bool makeNiFile(const std::string& niFile);

    /// @brief Point at a prebuilt NI index to load instead of re-reading the FASTA
    /// Purely an accelerator: the FASTA passed to the constructor stays authoritative,
    /// and an index that fails verification is dropped in favour of it.
    /// @param niFile NI file path supplied by the user
    void setNiFile(const std::string& niFile) { niIndexFile = niFile; }

    /// @brief Load reference genome compressed data from a user-specified NI file
    /// The payload is verified against the checksum recorded inside the file, so a
    /// damaged or tampered index is rejected rather than silently trusted.
    /// @param niFile NI file path supplied by the user
    /// @return true for success, false for failure
    bool loadFromNiFile(const std::string& niFile);

    /// @brief Tell whether a path holds a pbgz NI index, by content rather than by suffix
    /// @param file Path to probe
    /// @return true when the file carries the NI magic header
    static bool isNiFile(const std::string& file);

    /// @brief Build squash buffer in-memory directly from FASTA file
    /// No disk cache, no files written to user space.
    /// @return true for success, false for failure
    bool initSquashFromFasta();

    /// @brief Initialize reference genome compressed data through streaming
    /// Allocate buffer of specified size for storing compressed data
    /// @param squashLength Expected length of compressed data
    /// @return Pointer to allocated compressed data buffer
    uint8_t* initSquashByStream(int64_t squashLength);

    /// @brief Get ACTG base sequence of specified position range
    /// Extract ACTG base sequence of specified position and length from reference genome
    /// @param out Output buffer for storing ACTG characters
    /// @param outLength Length of output buffer
    /// @param actgPos Starting ACTG base position
    void getStretchActg(uint8_t* out, uint32_t outLength, uint64_t actgPos);

    /// @brief Get 2-bit encoding sequence of specified position range
    /// Extract 2-bit encoding sequence of specified position and length from reference genome (each character stores 2 bits)
    /// @param out Output buffer for storing 2-bit encoding
    /// @param outLength Length of output buffer
    /// @param actgPos Starting ACTG base position
    void getStretch2Bits1Char(uint8_t* out, uint32_t outLength, uint64_t actgPos);

    /// @brief Convert 2-bit encoding to ACTG character sequence
    /// Convert input 2-bit encoding data to corresponding ACTG characters
    /// @param src2Bit Input 2-bit encoding data
    /// @param src2BitsLen Length of input data
    /// @param dstActg Output buffer for storing ACTG characters
    void getActgFrom2Bits(const uint8_t* src2Bit, uint32_t src2BitsLen, uint8_t* dstActg);

    /// @brief Debug function: export hash table content to file
    /// Write hash table content to "hash_table" file for debugging and analysis
    void dumpHashTable();

    void dumpSquash();

    /// @brief Detect if file is in gzip format
    /// @param fileName File path to check
    /// @return true if gzip file, false otherwise
    bool isGzipFile(const std::string& fileName);

    /// @brief Create appropriate IOReader based on file type and hardware support
    /// @param fileName File path to read
    /// @return Unique pointer to IOReader instance
    std::unique_ptr<IOReader> createIOReader(const std::string& fileName);

private:
    /// @brief Check validity of reference genome parameters
    /// Verify key parameters like baseGroupStep and baseGroupLen meet requirements
    /// @return PBGZ_ERR_OK for success, error code for failure
    int32_t referencCheck();

    /// @brief Extract base groups from reference genome and calculate hash values
    /// Multi-threaded processing, split reference genome into base groups of baseGroupLen length with baseGroupStep step,
    /// calculate hash value for each base group and store in bgHash array
    /// @param bgHash Output parameter, pointer to array storing base group hash information
    void makeIndexFetchBaseGroup(BaseGroupHash *&bgHash);

    /// @brief Calculate hash table size and element count for each bucket
    /// Count elements in each hash bucket and allocate memory space for hash table
    /// @param hashTable Output parameter, vector containing hash table size information
    void makeIndexCalcHashTableSize(HashTable &hashTable);

    /// @brief Initialize hash table structure
    /// Initialize hash table memory layout based on calculated hash table size information,
    /// set starting position and conflict handling mechanism for each bucket
    /// @param hashTable Hash table size information
    /// @param hashBucketCurPos Output parameter, current position pointer array for each bucket
    void makeIndexInitHashTable(const HashTable& hashTable, uint32_t* &hashBucketCurPos);

    /// @brief Build hash table content
    /// Fill hash table with base group position information and handle hash conflicts
    /// @param bgHash Base group hash information array
    /// @param hashBucketCurPos Current position pointer array for each bucket
    void makeIndexBuildHashTable(const BaseGroupHash* bgHash, uint32_t* &hashBucketCurPos);

    /// @brief Sort hash table
    /// Sort position information in each bucket to optimize query performance
    void makeIndexSortHashTable();

private:
    /// @brief Calculate MD5 checksum of reference genome file
    /// @param refGeneFile Reference genome file path
    /// @param refGeneMd5 Output MD5 string
    /// @param refGeneLen Output file length
    void calculateReferenceMd5(const std::string& refGeneFile, std::string& refGeneMd5, int64_t& refGeneLen);

    /// @brief Write NI file header
    /// @param niWriter File writer for NI file
    /// @param offset Output offset for data section
    /// @return true for success, false for failure
    bool writeNiFileHeader(FileWriter& niWriter, int64_t& offset);

    /// @brief Process FASTA file and write compressed data
    /// @param refGeneFile Reference genome file path
    /// @param niWriter File writer for NI file
    /// @param md5 MD5 context for checksum calculation
    /// @param wlen Output written data length
    /// @return true for success, false for failure
    bool writeNiFile(const std::string& niFile);
    /// @brief Generate and write metadata to NI file
    /// @param niWriter File writer for NI file
    /// @param refGeneMd5 Reference genome MD5
    /// @param refGeneLen Reference genome length
    /// @param md5 MD5 context for data checksum
    /// @return true for success, false for failure
    bool writeNiFileMetadata(FileWriter& niWriter, const std::string& refGeneMd5, int64_t refGeneLen, std::string& md5);

private:
    // Progress bar
    GuardBar* guardBar;
    int64_t gbCurrent;
    int64_t gbTotal;

    // FASTA file checksum, currently using MD5
    std::string fastaChecksum;
    // FASTA file content length
    int64_t fastaLength;
    // Length of bases used as key when building index table, must be odd
    const uint32_t baseGroupLen = 31;
    // Base step when building index table
    const uint32_t baseGroupStep = 32;
    // Concurrency level
    uint32_t parallel;

    // Hash table information
    // Element count for each hash bucket
    uint32_t* hashBucketCnt;

    // Hash table buffer recording reference genome base position information
    uint32_t* hashTableBuffer;

    // Reference genome file
    std::string refGeneFile;

    // Optional prebuilt index; empty means always parse refGeneFile
    std::string niIndexFile;

    // Buffer for reference gene ACTG encoding
    uint8_t* refGeneSquash;

    // Corresponding length of reference gene ACTG encoding buffer
    int64_t refGeneSquashlen;

    // Buffer for matched bases after reference genome ACTG encoding
    uint8_t* refGeneSquashMatched;

    // Corresponding length
    uint64_t refGeneSquashMatchedlen;

    // ACTG hash information
    const int32_t hashBuckets = (32 << 20); // 32M

    // Mask corresponding to hash buckets
    const int32_t hashMask = hashBuckets - 1;

    // 1 byte squashed actg stretch to 4 actg
    const uint32_t actgStretch[256] = {
        0x41414141, 0x43414141, 0x54414141, 0x47414141, 0x41434141, 0x43434141, 0x54434141, 0x47434141, 0x41544141, 0x43544141, 0x54544141, 0x47544141, 0x41474141, 0x43474141, 0x54474141, 0x47474141,
        0x41414341, 0x43414341, 0x54414341, 0x47414341, 0x41434341, 0x43434341, 0x54434341, 0x47434341, 0x41544341, 0x43544341, 0x54544341, 0x47544341, 0x41474341, 0x43474341, 0x54474341, 0x47474341,
        0x41415441, 0x43415441, 0x54415441, 0x47415441, 0x41435441, 0x43435441, 0x54435441, 0x47435441, 0x41545441, 0x43545441, 0x54545441, 0x47545441, 0x41475441, 0x43475441, 0x54475441, 0x47475441,
        0x41414741, 0x43414741, 0x54414741, 0x47414741, 0x41434741, 0x43434741, 0x54434741, 0x47434741, 0x41544741, 0x43544741, 0x54544741, 0x47544741, 0x41474741, 0x43474741, 0x54474741, 0x47474741,
        0x41414143, 0x43414143, 0x54414143, 0x47414143, 0x41434143, 0x43434143, 0x54434143, 0x47434143, 0x41544143, 0x43544143, 0x54544143, 0x47544143, 0x41474143, 0x43474143, 0x54474143, 0x47474143,
        0x41414343, 0x43414343, 0x54414343, 0x47414343, 0x41434343, 0x43434343, 0x54434343, 0x47434343, 0x41544343, 0x43544343, 0x54544343, 0x47544343, 0x41474343, 0x43474343, 0x54474343, 0x47474343,
        0x41415443, 0x43415443, 0x54415443, 0x47415443, 0x41435443, 0x43435443, 0x54435443, 0x47435443, 0x41545443, 0x43545443, 0x54545443, 0x47545443, 0x41475443, 0x43475443, 0x54475443, 0x47475443,
        0x41414743, 0x43414743, 0x54414743, 0x47414743, 0x41434743, 0x43434743, 0x54434743, 0x47434743, 0x41544743, 0x43544743, 0x54544743, 0x47544743, 0x41474743, 0x43474743, 0x54474743, 0x47474743,
        0x41414154, 0x43414154, 0x54414154, 0x47414154, 0x41434154, 0x43434154, 0x54434154, 0x47434154, 0x41544154, 0x43544154, 0x54544154, 0x47544154, 0x41474154, 0x43474154, 0x54474154, 0x47474154,
        0x41414354, 0x43414354, 0x54414354, 0x47414354, 0x41434354, 0x43434354, 0x54434354, 0x47434354, 0x41544354, 0x43544354, 0x54544354, 0x47544354, 0x41474354, 0x43474354, 0x54474354, 0x47474354,
        0x41415454, 0x43415454, 0x54415454, 0x47415454, 0x41435454, 0x43435454, 0x54435454, 0x47435454, 0x41545454, 0x43545454, 0x54545454, 0x47545454, 0x41475454, 0x43475454, 0x54475454, 0x47475454,
        0x41414754, 0x43414754, 0x54414754, 0x47414754, 0x41434754, 0x43434754, 0x54434754, 0x47434754, 0x41544754, 0x43544754, 0x54544754, 0x47544754, 0x41474754, 0x43474754, 0x54474754, 0x47474754,
        0x41414147, 0x43414147, 0x54414147, 0x47414147, 0x41434147, 0x43434147, 0x54434147, 0x47434147, 0x41544147, 0x43544147, 0x54544147, 0x47544147, 0x41474147, 0x43474147, 0x54474147, 0x47474147,
        0x41414347, 0x43414347, 0x54414347, 0x47414347, 0x41434347, 0x43434347, 0x54434347, 0x47434347, 0x41544347, 0x43544347, 0x54544347, 0x47544347, 0x41474347, 0x43474347, 0x54474347, 0x47474347,
        0x41415447, 0x43415447, 0x54415447, 0x47415447, 0x41435447, 0x43435447, 0x54435447, 0x47435447, 0x41545447, 0x43545447, 0x54545447, 0x47545447, 0x41475447, 0x43475447, 0x54475447, 0x47475447,
        0x41414747, 0x43414747, 0x54414747, 0x47414747, 0x41434747, 0x43434747, 0x54434747, 0x47434747, 0x41544747, 0x43544747, 0x54544747, 0x47544747, 0x41474747, 0x43474747, 0x54474747, 0x47474747};

    // 1 byte squashed actg stretch to 4 2bits
    const uint32_t actgStretch2bits[256] = {
        0x00000000, 0x01000000, 0x02000000, 0x03000000, 0x00010000, 0x01010000, 0x02010000, 0x03010000, 0x00020000, 0x01020000, 0x02020000, 0x03020000, 0x00030000, 0x01030000, 0x02030000, 0x03030000,
        0x00000100, 0x01000100, 0x02000100, 0x03000100, 0x00010100, 0x01010100, 0x02010100, 0x03010100, 0x00020100, 0x01020100, 0x02020100, 0x03020100, 0x00030100, 0x01030100, 0x02030100, 0x03030100,
        0x00000200, 0x01000200, 0x02000200, 0x03000200, 0x00010200, 0x01010200, 0x02010200, 0x03010200, 0x00020200, 0x01020200, 0x02020200, 0x03020200, 0x00030200, 0x01030200, 0x02030200, 0x03030200,
        0x00000300, 0x01000300, 0x02000300, 0x03000300, 0x00010300, 0x01010300, 0x02010300, 0x03010300, 0x00020300, 0x01020300, 0x02020300, 0x03020300, 0x00030300, 0x01030300, 0x02030300, 0x03030300,
        0x00000001, 0x01000001, 0x02000001, 0x03000001, 0x00010001, 0x01010001, 0x02010001, 0x03010001, 0x00020001, 0x01020001, 0x02020001, 0x03020001, 0x00030001, 0x01030001, 0x02030001, 0x03030001,
        0x00000101, 0x01000101, 0x02000101, 0x03000101, 0x00010101, 0x01010101, 0x02010101, 0x03010101, 0x00020101, 0x01020101, 0x02020101, 0x03020101, 0x00030101, 0x01030101, 0x02030101, 0x03030101,
        0x00000201, 0x01000201, 0x02000201, 0x03000201, 0x00010201, 0x01010201, 0x02010201, 0x03010201, 0x00020201, 0x01020201, 0x02020201, 0x03020201, 0x00030201, 0x01030201, 0x02030201, 0x03030201,
        0x00000301, 0x01000301, 0x02000301, 0x03000301, 0x00010301, 0x01010301, 0x02010301, 0x03010301, 0x00020301, 0x01020301, 0x02020301, 0x03020301, 0x00030301, 0x01030301, 0x02030301, 0x03030301,
        0x00000002, 0x01000002, 0x02000002, 0x03000002, 0x00010002, 0x01010002, 0x02010002, 0x03010002, 0x00020002, 0x01020002, 0x02020002, 0x03020002, 0x00030002, 0x01030002, 0x02030002, 0x03030002,
        0x00000102, 0x01000102, 0x02000102, 0x03000102, 0x00010102, 0x01010102, 0x02010102, 0x03010102, 0x00020102, 0x01020102, 0x02020102, 0x03020102, 0x00030102, 0x01030102, 0x02030102, 0x03030102,
        0x00000202, 0x01000202, 0x02000202, 0x03000202, 0x00010202, 0x01010202, 0x02010202, 0x03010202, 0x00020202, 0x01020202, 0x02020202, 0x03020202, 0x00030202, 0x01030202, 0x02030202, 0x03030202,
        0x00000302, 0x01000302, 0x02000302, 0x03000302, 0x00010302, 0x01010302, 0x02010302, 0x03010302, 0x00020302, 0x01020302, 0x02020302, 0x03020302, 0x00030302, 0x01030302, 0x02030302, 0x03030302,
        0x00000003, 0x01000003, 0x02000003, 0x03000003, 0x00010003, 0x01010003, 0x02010003, 0x03010003, 0x00020003, 0x01020003, 0x02020003, 0x03020003, 0x00030003, 0x01030003, 0x02030003, 0x03030003,
        0x00000103, 0x01000103, 0x02000103, 0x03000103, 0x00010103, 0x01010103, 0x02010103, 0x03010103, 0x00020103, 0x01020103, 0x02020103, 0x03020103, 0x00030103, 0x01030103, 0x02030103, 0x03030103,
        0x00000203, 0x01000203, 0x02000203, 0x03000203, 0x00010203, 0x01010203, 0x02010203, 0x03010203, 0x00020203, 0x01020203, 0x02020203, 0x03020203, 0x00030203, 0x01030203, 0x02030203, 0x03030203,
        0x00000303, 0x01000303, 0x02000303, 0x03000303, 0x00010303, 0x01010303, 0x02010303, 0x03010303, 0x00020303, 0x01020303, 0x02020303, 0x03020303, 0x00030303, 0x01030303, 0x02030303, 0x03030303};

};
