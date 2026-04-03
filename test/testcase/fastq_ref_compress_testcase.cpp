/*
 * fastq_ref_compress_testcase.cpp - Test cases for FASTQ reference compression
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
#include <string>
#include <vector>
#include <cstdint>
#include <memory>

#define private public
#include "fastq_actuator.h"
#include <reference.h>
#include <io_wrapper.h>
#include <block_wrapper.h>
#undef private

namespace FastqRefTestData {
    const std::string testRefFile = "test_ref.fa";
    const uint32_t MAX_BLOCK_SIZE = 8 << 20;
};

class FastqRefCompressTest : public ::testing::Test {
public:
    void SetUp() override {
        generateReferenceFile(FastqRefTestData::testRefFile);

        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
        
        pInBlock = new RoughIOBlock(FastqRefTestData::MAX_BLOCK_SIZE);
        pInBlock->setMaxLineLen(80);
        pOutBlock = new RoughIOBlock(FastqRefTestData::MAX_BLOCK_SIZE);
        pReference = new Reference(FastqRefTestData::testRefFile, 1);
        pReference->makeIndex();
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
        if (pReference != nullptr) {
            delete pReference;
            pReference = nullptr;
        }

        std::remove(FastqRefTestData::testRefFile.c_str());
    }

    void generateReferenceFile(const std::string& filename) {
        std::ofstream fastaFile(filename);
        if (!fastaFile.is_open()) {
            return;
        }
        
        // Generate reference sequences with 80 characters per line
        
        // Chromosome 1: ATCG pattern repeated 20 times = 80 characters
        fastaFile << ">chr1\n";
        fastaFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
        fastaFile << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n";
        fastaFile << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";
        fastaFile << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";
        fastaFile << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n";
        fastaFile << "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n";
        
        // Chromosome 2: Sequence with N characters and mixed content
        fastaFile << ">chr2\n";
        fastaFile << "ATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGN\n";
        fastaFile << "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";
        fastaFile << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA\n";
        fastaFile << "ATNNCATNNCATNNCATNNCATNNCATNNCATNNCATNNCATNNCATNNCATNNCATNNCATNNCATNNCATNNCATNNC\n";
        
        // Chromosome 3: Short sequence
        fastaFile << ">chr3\n";
        fastaFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
        
        // Chromosome 4: Sequence with repetitive patterns
        fastaFile << ">chr4\n";
        fastaFile << "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT\n";
        fastaFile << "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC\n";
        fastaFile << "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA\n";
        
        // Chromosome 5: Sequence with edge cases (all Ns, mixed case)
        fastaFile << ">chr5\n";
        fastaFile << "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";
        fastaFile << "atcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcg\n";
        fastaFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
        
        fastaFile.close();
    }

protected:
    RoughIOBlock* pInBlock;
    RoughIOBlock* pOutBlock;
    Reference* pReference;
};

// Test basic functionality of mappingFastqGen2 function
TEST_F(FastqRefCompressTest, testMappingFastqGen2) {
    // Create FastqActuator instance
    FastqActuator actuator(pInBlock, pOutBlock, pReference);
    
    // Initialize encoder
    int32_t result = actuator.initEncoder();
    EXPECT_EQ(result, 0);

    uint32_t baseLength = 80;
    // Prepare output buffer
    uint8_t* outBuffer = new uint8_t[baseLength * 2];
    uint8_t* out = outBuffer;
    uint32_t outLength = 0;
    uint64_t mappingPos = 0;
    uint8_t mappingDir = 0;
    
    {
        // Prepare test data - a simple DNA sequence
        const uint8_t* base = (const uint8_t*)"ATCGATCGATCGATCGATCG";
        baseLength = strlen((char*)base);
        // const uint8_t expOut[] = {0x27, 0x27, 0x27, 0x27, 0x27};
        // Call mappingFastqGen2 function
        actuator.mappingFastqGen2(base, 20, out, outLength, mappingPos, mappingDir);
        
        // Verify result
        EXPECT_EQ(mappingPos, 0);
        EXPECT_EQ(mappingDir, 2);
        EXPECT_EQ(outLength, 20);
        // EXPECT_EQ(memcmp(out, expOut, outLength), 0);
    } 

    {
        const uint8_t* base = (const uint8_t*)"ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG";
        baseLength = strlen((char*)base);
         // Call mappingFastqGen2 function
        actuator.mappingFastqGen2(base, baseLength, out, outLength, mappingPos, mappingDir);

        fprintf(stderr, "outLength=%d, mappingPos=%ld, mappingDir=%d", outLength, mappingPos, mappingDir);
        EXPECT_EQ(mappingDir, 0);
    }
    
    // Clean up resources
    delete[] outBuffer;
}
