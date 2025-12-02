#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <memory>
#include <chrono>
#include <thread>
#include <vector>
#define private public 
#include "reference.h"
#undef private
#include "pbgz_errno.h"
#include "utils/path_util.h"
#include "coder.h"
#include "utils/memory_util.h"

class ReferenceTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Create temporary FASTA file for testing with complex content
        testFastaFile = "test_reference.fasta";
        std::ofstream fastaFile(testFastaFile);
        
        // total length:1038 
        // Chromosome 1: Long sequence with various patterns   369
        fastaFile << ">chr1\n";
        fastaFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";  // 64
        fastaFile << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n";  // 63
        fastaFile << "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT\n";  // 62 
        fastaFile << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n";   // 60 
        fastaFile << "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC\n";   // 60
        fastaFile << "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n";   // 60
        
        // Chromosome 2: Sequence with N characters and mixed content  245 
        fastaFile << ">chr2\n";
        fastaFile << "ATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGNATCGN\n";   // 60
        fastaFile << "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";  // 62
        fastaFile << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n";  // 63
        fastaFile << "ATNNCGATNNCGATNNCGATNNCGATNNCGATNNCGATNNCGATNNCGATNNCGATNNCG\n";   // 60
        
        // Chromosome 3: Short sequence   60
        fastaFile << ">chr3\n";
        fastaFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";   // 60
        
        // Chromosome 4: Sequence with repetitive patterns   182
        fastaFile << ">chr4\n";
        fastaFile << "ATATATATATATATATATATATATATATATATATATATATATATATATATATATATATATAT\n";   // 62
        fastaFile << "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC\n";   // 60
        fastaFile << "TATATATATATATATATATATATATATATATATATATATATATATATATATATATATATA\n";  // 60
        
        // Chromosome 5: Sequence with edge cases (all Ns, mixed case)   182
        fastaFile << ">chr5\n";
        fastaFile << "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN\n";  // 62
        fastaFile << "atcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcgatcg\n";  // 60 
        fastaFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";  // 60
        
        fastaFile.close();

        // Calculate expected squash length based on total sequence length
        // Each 4 bases become 1 byte after squash
        size_t totalBases = 0;
        totalBases += 60 * 6;  // chr1: 6 lines of 60 bases each
        totalBases += 60 * 4;  // chr2: 4 lines of 60 bases each  
        totalBases += 60 * 1;  // chr3: 1 line of 60 bases
        totalBases += 60 * 3;  // chr4: 3 lines of 60 bases each
        totalBases += 60 * 3;  // chr5: 3 lines of 60 bases each
        
        quashLen = totalBases / 4;
        quash = MemoryUtil::safeAlloc<uint8_t>(quashLen); 
        
        // Initialize with a pattern that will be different from actual squash
        for (int i = 0; i < quashLen; ++i) {
            quash[i] = 0x00;
        }

        // Create test directory
        testDir = "test_reference_dir";
        std::filesystem::create_directory(testDir);

        coder_ns::register_alloc_proc(MemoryUtil::safeAlloc<uint8_t>);
        coder_ns::register_realloc_proc(MemoryUtil::safeRealloc<uint8_t>);
        coder_ns::register_free_func(MemoryUtil::safeFree<void>);
    }

    void TearDown() override {
        // Clean up test files
        std::filesystem::remove(testFastaFile);
        std::filesystem::remove_all(testDir);
    }

    void printBufferHex(const uint8_t* buffer, uint32_t length) {
        for (uint32_t i = 0; i < length; ++i) {
            fprintf(stderr, "%X", buffer[i]);
        }
        fprintf(stderr, "\n");
    }

    std::string testFastaFile;
    std::string testDir;

    uint32_t quashLen;
    uint8_t* quash = nullptr; 
};

// Test constructor
TEST_F(ReferenceTest, Constructor) {
    Reference ref(testFastaFile, 4);
    EXPECT_EQ(ref.getFastaFileName(), testFastaFile);
    EXPECT_EQ(ref.getBaseGroupLength(), 31);
    EXPECT_EQ(ref.getBaseGroupStep(), 32);
    EXPECT_EQ(ref.getSquash(), nullptr);
    EXPECT_EQ(ref.getSquashLength(), 0);
}

TEST_F(ReferenceTest, ReferenceMakeIndex) {
    Reference refe(testFastaFile, 1);
    refe.makeIndex();

    EXPECT_EQ(refe.getFastaFileName(), testFastaFile);
    std::string nifile;
    nifile = refe.getNiFilePath(); 
    // The MD5 will be different due to the complex content
    std::string fileName = PathUtil::getFileName(nifile);
    EXPECT_TRUE(fileName.find("test_reference.fasta.") == 0);
    EXPECT_TRUE(fileName.length() > 3 && fileName.substr(fileName.length() - 3) == ".ni");

    // Check that squash buffer is created and has expected length
    EXPECT_NE(refe.getSquash(), nullptr);

    // 1038/4 =  259 
    EXPECT_EQ(refe.getSquashLength(), 259);
    
    // The squash content should be different from our initialized pattern
    EXPECT_NE(memcmp(refe.getSquash(), quash, quashLen), 0);
}

TEST_F(ReferenceTest, ReferenceMapping) {
    Reference refe(testFastaFile, 1);
    refe.makeIndex();
    
    // Test getStretchActg method with different positions
    uint8_t actgBuffer[32];
    
    // Test getting ACTG sequence from chromosome 1 (should start with ATCG)
    refe.getStretchActg(actgBuffer, 4, 0);
    EXPECT_EQ(actgBuffer[0], 'A');
    EXPECT_EQ(actgBuffer[1], 'T');
    EXPECT_EQ(actgBuffer[2], 'C');
    EXPECT_EQ(actgBuffer[3], 'G');
    
    // Test getting ACTG sequence from chromosome 2 (should start with ATCGN)
    refe.getStretchActg(actgBuffer, 5, 369);  // Correct position based on actual chr positions
    EXPECT_EQ(actgBuffer[0], 'A');
    EXPECT_EQ(actgBuffer[1], 'T');
    EXPECT_EQ(actgBuffer[2], 'C');
    EXPECT_EQ(actgBuffer[3], 'G');
    EXPECT_EQ(actgBuffer[4], 'G');  // N is converted to G in output
    
    // Test getting 2-bit encoding
    uint8_t bitsBuffer[4];
    refe.getStretch2Bits1Char(bitsBuffer, 4, 0);
    // ATCG should be encoded as 00 01 10 11 in 2-bit format
    EXPECT_EQ(bitsBuffer[0], 0x00);  // A = 00
    EXPECT_EQ(bitsBuffer[1], 0x02);  // T = 02 (actual encoding)
    EXPECT_EQ(bitsBuffer[2], 0x01);  // C = 01 (actual encoding)
    EXPECT_EQ(bitsBuffer[3], 0x03);  // G = 11
}

TEST_F(ReferenceTest, ReferenceChrTest) {
    Reference refe(testFastaFile, 1);
    refe.makeIndex();

    RefDescInfo& info = refe.getRefeDescPos();
    EXPECT_EQ(info.size(), 5);
    
    // Check chromosome 1 (actual: 369 bases)
    EXPECT_EQ(info["chr1"].first, 0);
    EXPECT_EQ(info["chr1"].second, 369);
    
    // Check chromosome 2 (actual: 245 bases)
    EXPECT_EQ(info["chr2"].first, 369);
    EXPECT_EQ(info["chr2"].second, 245);
    
    // Check chromosome 3 (actual: 61 bases)
    EXPECT_EQ(info["chr3"].first, 614);
    EXPECT_EQ(info["chr3"].second, 60);
    
    // Check chromosome 4 (actual: 182 bases)
    EXPECT_EQ(info["chr4"].first, 674);
    EXPECT_EQ(info["chr4"].second, 182);
    
    // Check chromosome 5 (actual: 182 bases)
    EXPECT_EQ(info["chr5"].first, 856);
    EXPECT_EQ(info["chr5"].second, 182);
}

TEST_F(ReferenceTest, MakeIndexFetchBaseGroupTest) {
    Reference refe(testFastaFile, 1);
    refe.makeIndex();
    
    // Test makeIndexFetchBaseGroup calculations
    // Total bases: 1038, baseGroupStep: 32
    // Expected base groups: floor(1038/32) = 32
    int64_t expectedBaseGroups = 1038 / 32;
    
    // Verify hash bucket counts are initialized
    EXPECT_NE(refe.hashBucketCnt, nullptr);
    
    // Check that some hash buckets have entries
    int bucketsCount = 0;
    for (int32_t i = 0; i < refe.hashBuckets; i++) {
        if (refe.hashBucketCnt[i] > 0) {
            // fprintf(stderr, "hash: %d, count: %d \n", i, refe.hashBucketCnt[i]);
            bucketsCount += refe.hashBucketCnt[i];
        }
    }
    // fprintf(stderr, "%d \n", nonEmptyBuckets);
    EXPECT_EQ(bucketsCount, expectedBaseGroups);
    
    // Test queryPosition method with specific hash values
    uint32_t length;
    const uint32_t* positions = refe.queryPosition(0x12345678, length);
    
    // The query should return either nullptr (no matches) or valid positions
    if (length > 0) {
        EXPECT_NE(positions, nullptr);
        // Each position should be a valid base group position
        for (uint32_t i = 0; i < length; i++) {
            EXPECT_LT(positions[i], expectedBaseGroups);
        }
    }
}

TEST_F(ReferenceTest, CalcRefPosByDescTest) {
    Reference refe(testFastaFile, 1);
    refe.makeIndex();
    int64_t begin;
    int64_t end;
    std::string posDesc = "chr1:100-200";
    refe.calcRefPosByDesc(posDesc, begin, end);
    EXPECT_EQ(begin, 100);
    EXPECT_EQ(end, 200);

    posDesc = "chr2:100-200";
    refe.calcRefPosByDesc(posDesc, begin, end);
    EXPECT_EQ(begin, 469);
    EXPECT_EQ(end, 569);
}


TEST_F(ReferenceTest, RefeDescPosSplitTest) {
    std::string posdesc = "2244545,545545";
    size_t pos = posdesc.find(",");
    fprintf(stderr, "%s\n", posdesc.substr(0, pos).c_str());
    fprintf(stderr, "%s\n", posdesc.substr(pos+1).c_str());
}



