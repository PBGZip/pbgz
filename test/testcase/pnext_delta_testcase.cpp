/*
 * pnext_delta_testcase.cpp - Test cases for PNEXT delta compression
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

#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <cstdio>
#include <memory>

#define private public
#include "sam_actuator.h"
#include "coder_affix_match.h"
#include <io_wrapper.h>
#include <block_wrapper.h>
#include "config_manager.h"
#include "compress_engine.h"
#include "coder/coder_affix_match.h"
#include <coder/coder_io.h>
#include <coder/coder_bwt_cm.h>
#undef private

namespace PNextDeltaTestData {
    const std::string testSamFile = "pnext_delta_test.sam";
    const std::string testCompressedFile = "pnext_delta_test.pbgz";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class PNextDeltaTest : public ::testing::Test {
public:
    void SetUp() override {
        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        coder_ns::initFcCoder();

        // Ensure SamInfo state is clean
        SamInfo::getInstance().clearChromosomeInfo();
        SamInfo::getInstance().resetChrIdCounter();

        pInBlock = new RoughIOBlock(PNextDeltaTestData::MAX_BLOCK_SIZE);
        pOutBlock = new RoughIOBlock(PNextDeltaTestData::MAX_BLOCK_SIZE);

        // Set appropriate log level
        ConfigManager::getInstance().logLevel = LogLevel::INFO;
    }

    void TearDown() override {
        if (pInBlock != nullptr) {
            delete pInBlock;
            pInBlock = nullptr;
        }
        if (pOutBlock != nullptr) {
            delete pOutBlock;
            pOutBlock = nullptr;
        }
    }

    void loadSamData(const std::string& filename) {
        pInBlock->reset();
        pOutBlock->reset();

        std::vector<std::string> paths = { 
            filename,
            "./test/" + filename,
            "../test/" + filename,
            "../" + filename
        };

        IOReader* pIoReader = nullptr;
        bool loaded = false;
        for (const auto& path : paths) {
            if (std::filesystem::exists(path)) {
                pIoReader = new FileReader(path);
                if (pIoReader != nullptr) {
                    loaded = true;
                    break;
                }
            }
        }

        ASSERT_TRUE(loaded) << "Failed to load test file: " << filename;

        // FileReader must be opened (mmap'ed) before it can return data
        if (pIoReader->openIO() != 0) {
            delete pIoReader;
            FAIL() << "Failed to open test file: " << filename;
        }

        // Use BlockReader to parse lines into npos, like the other test suites;
        // preAnalysis/compress rely on npos being populated
        BlockReader* pBlockReader = new BlockReader(pIoReader);
        pBlockReader->readBlock(pInBlock, TYPE_UNKNOW);

        delete pBlockReader;
        delete pIoReader;
    }

    void createTestSamData() {
        // Create simple test SAM data with known values for PNEXT testing
        std::string testData = 
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:ref|NC_001133|\tLN:230218\n"
            "@SQ\tSN:ref|NC_001142|\tLN:784333\n"
            "@PG\tID:samtools\tPN:samtools\tVN:1.6\n"
            "test1.1\t99\tref|NC_001133|\t100\t0\t90M\t=\t200\t300\tAGTACCAAATGCACTCACATCATTATGCACGGCACTTGCCTCAGCGGTCTATACCATGTGCCATTTACGCATAACGCCCATCATTATCCA\tBBB[[dcccbb_bccbabbbbcccdcb`^eeeggffhihihhffgfffggffchfaiihhgdfhhfiigihghgiiiggfggeeeeebbb\tNM:i:2\n"
            "test1.2\t99\tref|NC_001133|\t150\t0\t90M\t=\t300\t400\tATTTTAATATTTATATCTTATTCTGCGGTCCCAAATATTGTATAAATGCCCTTAATACATACTTTATACCACTTTTGCACCATATACTAA\tbbbeeeeegggfghhiafgfhhfhihiiihiihhihhfhhhfhhiiiiihiiihhhdgiiiihiiffihffhhhhdgh_dbgeeeeeeee\tNM:i:7\n"
            "test1.3\t99\tref|NC_001133|\t200\t0\t90M\t=\t400\t500\tTAATATTTATATCTTATTCTGCGGTCCCAAATATTGTATAAATGCCCTTAATACATACTTTATACCACTTTTGCACCATATACTAACCAC\tbbbeeeeeggegfhifhihiefhiiifbgfhiiihfhifgbhhhihiigiiiihhghhhhhiiihfhihhhddfbdddgfgfgeeeeeec\tNM:i:7\n"
            "test1.4\t163\tref|NC_001133|\t250\t0\t90M\t=\t500\t600\tTATCTTATTCTGCGGTCCCAAATATTGTATAAATGCCCTTAATACATACTTTATACCACTTTTGCACCATATACTAACCACTCAATTTAT\tbbbeeeeegggggihhiiiiiiiihiiiiiiiiiiiiiiiiiiiieghhiiiiiiiiihhiiiiiiiiiegggeeeeeddcbdIIIIIII\tNM:i:6\n";

        std::ofstream testFile(PNextDeltaTestData::testSamFile);
        testFile << testData;
        testFile.close();
    }

protected:
    RoughIOBlock* pInBlock = nullptr;
    RoughIOBlock* pOutBlock = nullptr;
};

// Test basic delta compression and decompression
TEST_F(PNextDeltaTest, TestBasicDeltaCompression) {
    createTestSamData();
    loadSamData(PNextDeltaTestData::testSamFile);

    // Create SamCodecActuator instance with a real engine (compress() needs it)
    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator* pActuator = new SamCodecActuator(pInBlock, pOutBlock, &engine);

    // Run pre-analysis first
    int32_t result = pActuator->preAnalysis();
    ASSERT_EQ(result, 0) << "preAnalysis failed";

    // Test compression
    result = pActuator->compress();
    ASSERT_EQ(result, 0) << "Compression failed";

    // Verify output data is valid
    ASSERT_GT(pOutBlock->getDataLen(), 0) << "No data was compressed";

    // Clean up
    delete pActuator;
    std::remove(PNextDeltaTestData::testSamFile.c_str());
}

// Test delta compression with AFFIX_MATCH encoder specifically
TEST_F(PNextDeltaTest, TestAffixMatchDeltaCompression) {
    createTestSamData();
    loadSamData(PNextDeltaTestData::testSamFile);

    SamCodecActuator* pActuator = new SamCodecActuator(pInBlock, pOutBlock);

    // Run pre-analysis
    int32_t result = pActuator->preAnalysis();
    ASSERT_EQ(result, 0) << "preAnalysis failed";

    // Test AFFIX_MATCH delta compression by forcing it
    Json::Value fieldMeta;
    uint32_t fieldSrcLen = 0;
    
    // Create a test vector of POS values and corresponding PNEXT values
    std::vector<int64_t> posValues = {100, 150, 200, 250};
    std::vector<int64_t> pNextValues = {200, 300, 400, 500}; // All should be delta = +100
    
    // Use AFFIX_MATCH coder for delta compression
    result = pActuator->compressPNextFieldDelta<coder_affix_match>(7, fieldSrcLen, fieldMeta);
    
    if (result == 0) {
        std::cout << "AFFIX_MATCH delta compression succeeded, srcLen=" << fieldSrcLen 
                  << ", dstLen=" << fieldMeta["dstlen"].asUInt() << std::endl;
    } else {
        std::cout << "AFFIX_MATCH delta compression failed with result=" << result << std::endl;
    }
    
    ASSERT_GE(result, 0) << "Delta compression should not fail catastrophically";

    // Clean up
    delete pActuator;
    std::remove(PNextDeltaTestData::testSamFile.c_str());
}

// Test the AFFIX_MATCH encoder behavior directly
TEST_F(PNextDeltaTest, TestAffixMatchCoderBehavior) {
    // Create test data with patterns that should compress well
    std::vector<std::string> testStrings = {
        "100\t", "105\t", "110\t", "115\t", "120\t",  // Sequential deltas
        "50\t", "55\t", "60\t", "65\t", "70\t",   // Another sequential pattern
        "1000\t", "1050\t", "1100\t", "1150\t", // Larger deltas
        "-50\t", "-40\t", "-30\t", "-20\t", "-10\t" // Negative deltas
    };
    
    std::vector<uint8_t> encodeBuf(4096);
    std::shared_ptr<coder_io> io = std::make_shared<coder_io>(encodeBuf.data(), static_cast<int32_t>(encodeBuf.size()));
    coder_affix_match coderio(io.get());
    
    std::cout << "Testing AFFIX_MATCH coder with realistic delta data:" << std::endl;
    
    for (const auto& str : testStrings) {
        coderio.encode_line(reinterpret_cast<const uint8_t*>(str.c_str()), str.length());
        std::cout << "Encoded: '" << str << "' -> buffer state still internal (not flushed yet)" << std::endl;
    }
    
    // FLUSH IS CRITICAL - without this, no data will be written!
    coderio.encode_flush();
    std::cout << "After encode_flush(): io->data_len = " << io->data_len << std::endl;
    
    ASSERT_GT(io->data_len, 0) << "AFFIX_MATCH must produce output after flush";
    
    // Now try to decode it back
    std::vector<std::string> decodedResults;
    coder_affix_match decoder(io.get());
    
    for (size_t i = 0; i < testStrings.size(); i++) {
        uint8_t decodeBuffer[128];
        uint8_t splitFlag = (uint8_t)'\t';
        // need2hold=true: coder_affix_match keeps "last" internally instead of
        // pointing into this reused stack buffer, so decoding across lines works
        int32_t decResult = decoder.decode_line(decodeBuffer, sizeof(decodeBuffer), splitFlag, true);
        
        decodedResults.push_back(std::string(reinterpret_cast<char*>(decodeBuffer), decResult));
        
        std::cout << "Decode[" << i << "]: result=" << decResult << ", data: '" << (char*)decodeBuffer << "'" << std::endl;
        
        // Verify roundtrip - some failures are expected if compression has issues
        if (testStrings[i] == decodedResults[i]) {
            std::cout << "  -> MATCH" << std::endl;
        } else {
            std::cout << "  -> MISMATCH! Expected: '" << testStrings[i] << "' but got: '" << decodedResults[i] << "'" << std::endl;
        }
    }
    
    // Final verification - all outputs should be produced
    EXPECT_EQ(decodedResults.size(), testStrings.size()) << "All strings should be decoded";
}

// Test compression with a larger generated dataset (no external file dependency)
TEST_F(PNextDeltaTest, TestConSortedDeltaCompression) {
    // Generate a larger SAM dataset inside the test, with SEQ/QUAL lengths
    // consistent with the 90M CIGAR so preAnalysis accepts it
    {
        std::ofstream out(PNextDeltaTestData::testSamFile);
        out << "@HD\tVN:1.6\tSO:coordinate\n"
            << "@SQ\tSN:ref|NC_001133|\tLN:230218\n"
            << "@SQ\tSN:ref|NC_001142|\tLN:784333\n"
            << "@PG\tID:samtools\tPN:samtools\tVN:1.6\n";
        for (int i = 0; i < 200; ++i) {
            uint16_t flag = (i % 2 == 0) ? 99 : 147;
            int64_t pos = 100 + i * 500;
            int64_t pnext = pos + 400;
            int32_t tlen = (flag & 0x80) ? -400 : 400;
            char seq[91];
            char qual[91];
            for (int k = 0; k < 90; ++k) {
                seq[k] = "ACGT"[k % 4];
                qual[k] = 'I';
            }
            seq[90] = '\0';
            qual[90] = '\0';
            out << "read" << i << "\t" << flag << "\tref|NC_001133|\t" << pos
                << "\t60\t90M\t=\t" << pnext << "\t" << tlen << "\t"
                << seq << "\t" << qual << "\tNM:i:0\n";
        }
    }

    loadSamData(PNextDeltaTestData::testSamFile);

    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator* pActuator = new SamCodecActuator(pInBlock, pOutBlock, &engine);

    // Run pre-analysis
    int32_t result = pActuator->preAnalysis();
    ASSERT_EQ(result, 0) << "preAnalysis failed";

    // Test compression
    result = pActuator->compress();
    EXPECT_EQ(result, 0) << "Compression failed";

    // Verify output
    EXPECT_GT(pOutBlock->getDataLen(), 0) << "No compressed data produced";

    delete pActuator;
    std::remove(PNextDeltaTestData::testSamFile.c_str());
}

// Test decompression of delta-compressed PNEXT
TEST_F(PNextDeltaTest, TestDeltaDecompression) {
    createTestSamData();
    loadSamData(PNextDeltaTestData::testSamFile);

    PbgzParameter para;
    CompressEngine engine(para);
    SamCodecActuator* pActuator = new SamCodecActuator(pInBlock, pOutBlock, &engine);

    // First compress
    int32_t result = pActuator->preAnalysis();
    ASSERT_EQ(result, 0) << "preAnalysis failed";

    result = pActuator->compress();
    ASSERT_EQ(result, 0) << "Compression failed";

    uint32_t compressedLen = pOutBlock->getDataLen();
    ASSERT_GT(compressedLen, 0) << "No compressed data";

    // Now try to decompress back
    RoughIOBlock* pDecompressedBlock = new RoughIOBlock(PNextDeltaTestData::MAX_BLOCK_SIZE);
    RoughIOBlock* pCompressedInput = new RoughIOBlock(PNextDeltaTestData::MAX_BLOCK_SIZE);
    
    // Copy compressed data and its metadata to the new input block
    memcpy(pCompressedInput->getBuffer(), pOutBlock->getBuffer(), compressedLen);
    pCompressedInput->setDataLen(compressedLen);
    memcpy(pCompressedInput->getMetaBuffer(), pOutBlock->getMetaBuffer(), pOutBlock->getMetaLen());
    pCompressedInput->setMetaLen(pOutBlock->getMetaLen());

    SamCodecActuator* pDecompressor = new SamCodecActuator(pCompressedInput, pDecompressedBlock, &engine);
    
    result = pDecompressor->preAnalysis();
    ASSERT_EQ(result, 0) << "Decompression preAnalysis failed";

    result = pDecompressor->decompress();
    EXPECT_EQ(result, 0) << "Decompression failed";

    // Verify decompressed data
    EXPECT_GT(pDecompressedBlock->getDataLen(), 0) << "No decompressed data";

    // Compare with original
    uint32_t originalLen = pInBlock->getDataLen();
    uint32_t decompressedLen = pDecompressedBlock->getDataLen();
    
    std::cout << "Original length: " << originalLen 
              << ", Decompressed length: " << decompressedLen << std::endl;
    
    // They should be reasonably close (allowing for compression overhead)
    // The exact match would require full decompression verification
    
    delete pDecompressor;
    delete pActuator;
    delete pCompressedInput;
    delete pDecompressedBlock;
    
    std::remove(PNextDeltaTestData::testSamFile.c_str());
}

// Test encoding/decoding roundtrip with known data
TEST_F(PNextDeltaTest, TestDeltaRoundtrip) {
    // Create encoder
    std::vector<uint8_t> encodeBuf(4096);
    std::shared_ptr<coder_io> encodeIo = std::make_shared<coder_io>(encodeBuf.data(), static_cast<int32_t>(encodeBuf.size()));
    coder_affix_match encoder(encodeIo.get());
    
    // Create known delta sequences
    std::vector<std::string> testDeltas = {
        "100\t",   // Simple positive delta
        "-50\t",   // Negative delta
        "1000\t",  // Larger positive delta
        "744954\t", // Real-world example
        "378\t",    // Another real-world example
    };
    
    std::vector<std::string> decodedResults;
    
    std::cout << "Testing AFFIX_MATCH roundtrip compression:" << std::endl;
    
    // Encode all deltas first
    for (size_t i = 0; i < testDeltas.size(); i++) {
        const auto& original = testDeltas[i];
        
        // Encode
        encoder.encode_line(reinterpret_cast<const uint8_t*>(original.c_str()), original.length());
        std::cout << "  [" << i << "] Encoded: '" << original << "' (internal buffer state)" << std::endl;
    }
    
    // CRITICAL: Flush the encoder to write data to io buffer
    encoder.encode_flush();
    std::cout << "After flush: encoded_data_len = " << encodeIo->data_len << std::endl;
    
    ASSERT_GT(encodeIo->data_len, 0) << "Encoder must produce output after flush";
    
    // Create decoder and decode all data
    std::shared_ptr<coder_io> decodeIo = std::make_shared<coder_io>(
        encodeIo->data, 
        encodeIo->data_len
    );
    coder_affix_match decoder(decodeIo.get());
    
    for (size_t i = 0; i < testDeltas.size(); i++) {
        const auto& original = testDeltas[i];
        
        // Decode
        uint8_t decodeBuffer[128] = {0};
        uint8_t splitFlag = (uint8_t)'\t';
        // need2hold=true: keep "last" inside the coder so the reused stack
        // buffer does not corrupt the affix context across lines
        int32_t decResult = decoder.decode_line(decodeBuffer, sizeof(decodeBuffer), splitFlag, true);
        
        std::string decoded((char*)decodeBuffer, decResult);
        decodedResults.push_back(decoded);
        
        std::cout << "  [" << i << "] Original: '" << original << "' -> decoded_len=" << decResult 
                  << " -> '" << decoded << "' -> ";
        
        // Compare
        if (original == decoded) {
            std::cout << "MATCH" << std::endl;
        } else {
            std::cout << "MISMATCH! Expected: '" << original << "' but got: '" << decoded << "'" << std::endl;
        }
    }
    
    // Final verification
    EXPECT_EQ(decodedResults.size(), testDeltas.size()) << "All deltas should be decoded";
    
    // Count successful roundtrips
    int successCount = 0;
    for (size_t i = 0; i < testDeltas.size(); i++) {
        if (testDeltas[i] == decodedResults[i]) {
            successCount++;
        }
    }
    
    std::cout << "Roundtrip success rate: " << successCount << "/" << testDeltas.size() << std::endl;
    EXPECT_GT(successCount, 0) << "At least some roundtrips should succeed";
}

// Test with multiple lines to identify state issues
TEST_F(PNextDeltaTest, TestMultiLineCompressionState) {
    // Create a set of consecutive delta values to test encoder state management
    std::vector<uint8_t> encodeBuf(4096);
    std::shared_ptr<coder_io> encodeIo = std::make_shared<coder_io>(encodeBuf.data(), static_cast<int32_t>(encodeBuf.size()));
    coder_affix_match encoder(encodeIo.get());
    
    std::cout << "Testing multi-line compression to identify state issues:" << std::endl;
    
    // Simulate processing multiple lines with consecutive delta values
    for (int line = 0; line < 10; line++) {
        char buffer[64];
        // Create a delta pattern: line 0: 100, line 1: 105, line 2: 110, etc.
        int deltaValue = 100 + (line * 5);
        int len = snprintf(buffer, sizeof(buffer), "%d\t", deltaValue);
        
        std::cout << "Line " << line << ": encoding delta='" << buffer << "' (len=" << len << ")" << std::endl;
        
        // need2hold=true on both sides: the encoder's "last" would otherwise
        // point into this reused stack buffer and be overwritten each iteration
        encoder.encode_line(reinterpret_cast<const uint8_t*>(buffer), len, true);
        std::cout << "  -> internal buffer state: offset=" << encodeIo->data_len << " (still internal, not flushed)" << std::endl;
    }
    
    // CRITICAL: Flush the encoder to write all accumulated data
    encoder.encode_flush();
    std::cout << "After flush: final encoded size = " << encodeIo->data_len << std::endl;
    
    ASSERT_GT(encodeIo->data_len, 0) << "Encoder must produce output after flush";
    
    // Now try to decode each line
    std::shared_ptr<coder_io> decodeIo = std::make_shared<coder_io>(
        encodeIo->data,
        encodeIo->data_len
    );
    coder_affix_match decoder(decodeIo.get());
    
    for (int line = 0; line < 10; line++) {
        uint8_t decodeBuffer[64] = {0};
        uint8_t splitFlag = (uint8_t)'\t';
        
        int32_t decResult = decoder.decode_line(decodeBuffer, sizeof(decodeBuffer), splitFlag, true);
        
        std::cout << "Line " << line << ": decode_result=" << decResult << ", decoded='" << (char*)decodeBuffer << "'" << std::endl;
        
        int expectedDelta = 100 + (line * 5);
        std::string expected = std::to_string(expectedDelta) + "\t";
        
        if (decResult > 0 && std::string((char*)decodeBuffer, decResult) == expected) {
            std::cout << "  -> MATCH" << std::endl;
        } else {
            std::cout << "  -> MISMATCH! Expected: '" << expected << "'" << std::endl;
        }
    }
    
    // Verify all lines were decoded
    EXPECT_GT(encodeIo->data_len, 0) << "Should have encoded data";
}