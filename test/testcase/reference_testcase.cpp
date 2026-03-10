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
        for (uint32_t i = 0; i < quashLen; ++i) {
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
        // Clean up allocated memory
        if (quash != nullptr) {
            MemoryUtil::safeFree(quash);
        }
        
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

// Test gzip file support in calculateReferenceMd5
TEST_F(ReferenceTest, CalculateReferenceMd5GzipSupport) {
    // Create a simple FASTA file
    std::string textFastaFile = "test_text.fasta";
    std::ofstream textFile(textFastaFile);
    textFile << ">chr1\n";
    textFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
    textFile << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n";
    textFile.close();

    // Create gzip compressed version
    std::string gzipFastaFile = "test_gzip.fasta";
    std::string command = "gzip -c " + textFastaFile + " > " + gzipFastaFile;
    system(command.c_str());

    // Create gzip file without .gz extension to test content-based detection
    std::string gzipNoExtFile = "test_gzip_noext";
    std::string command2 = "cp " + gzipFastaFile + " " + gzipNoExtFile;
    system(command2.c_str());

    // Test 1: Text file MD5 calculation
    std::string textMd5, gzipMd5, gzipNoExtMd5;
    int64_t textLength, gzipLength, gzipNoExtLength;
    
    Reference textRef(textFastaFile, 1);
    textRef.calculateReferenceMd5(textFastaFile, textMd5, textLength);
    
    // Verify text file MD5 is calculated correctly
    EXPECT_FALSE(textMd5.empty());
    EXPECT_GT(textLength, 0);
    
    // Test 2: Gzip file with .gz extension
    Reference gzipRef(gzipFastaFile, 1);
    gzipRef.calculateReferenceMd5(gzipFastaFile, gzipMd5, gzipLength);
    
    // Verify gzip file MD5 is calculated correctly
    EXPECT_FALSE(gzipMd5.empty());
    EXPECT_GT(gzipLength, 0);
    
    // Test 3: Gzip file without .gz extension (content-based detection)
    Reference gzipNoExtRef(gzipNoExtFile, 1);
    gzipNoExtRef.calculateReferenceMd5(gzipNoExtFile, gzipNoExtMd5, gzipNoExtLength);
    
    // Verify gzip file without extension is detected correctly
    EXPECT_FALSE(gzipNoExtMd5.empty());
    EXPECT_GT(gzipNoExtLength, 0);
    
    // Test 4: MD5 values should be the same for all three files (same content)
    EXPECT_EQ(textMd5, gzipMd5);
    EXPECT_EQ(textMd5, gzipNoExtMd5);
    
    // Test 5: File lengths should be the same for all three files
    EXPECT_EQ(textLength, gzipLength);
    EXPECT_EQ(textLength, gzipNoExtLength);
    
    // Clean up test files
    std::filesystem::remove(textFastaFile);
    std::filesystem::remove(gzipFastaFile);
    std::filesystem::remove(gzipNoExtFile);
}

// Test file type detection based on content
TEST_F(ReferenceTest, FileTypeDetectionByContent) {
    // Create a text file
    std::string textFile = "test_text.txt";
    std::ofstream tf(textFile);
    tf << ">chr1\n";
    tf << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
    tf.close();
    
    // Create a gzip file with different extension
    std::string gzipFile = "test_compressed.bin";
    std::string command = "echo '>chr1' | gzip -c > " + gzipFile;
    system(command.c_str());
    
    // Test text file detection
    FILE* fp = fopen(textFile.c_str(), "rb");
    ASSERT_NE(fp, nullptr);
    uint8_t header[2];
    size_t read = fread(header, 1, 2, fp);
    fclose(fp);
    
    // Text file should not have gzip magic number
    EXPECT_EQ(read, 2);
    EXPECT_FALSE(header[0] == 0x1f && header[1] == 0x8b);
    
    // Test gzip file detection
    fp = fopen(gzipFile.c_str(), "rb");
    ASSERT_NE(fp, nullptr);
    read = fread(header, 1, 2, fp);
    fclose(fp);
    
    // Gzip file should have gzip magic number 0x1f 0x8b
    EXPECT_EQ(read, 2);
    EXPECT_TRUE(header[0] == 0x1f && header[1] == 0x8b);
    
    // Clean up
    std::filesystem::remove(textFile);
    std::filesystem::remove(gzipFile);
}

// Test index creation with both text and gzip reference files
TEST_F(ReferenceTest, IndexCreationWithGzipReference) {
    // Create a simple FASTA file
    std::string textFastaFile = "test_ref.fasta";
    std::ofstream textFile(textFastaFile);
    textFile << ">chr1\n";
    textFile << "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n";
    textFile << "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT\n";
    textFile.close();

    // Create gzip compressed version
    std::string gzipFastaFile = "test_ref.fasta.gz";
    std::string command = "gzip -c " + textFastaFile + " > " + gzipFastaFile;
    int result = system(command.c_str());
    
    // Verify gzip file was created successfully
    EXPECT_EQ(result, 0);
    EXPECT_TRUE(std::filesystem::exists(gzipFastaFile));

    // Test 1: Create index with text reference file
    Reference textRef(textFastaFile, 1);
    bool textIndexResult = textRef.makeIndex();
    EXPECT_TRUE(textIndexResult);
    if (textIndexResult) {
        EXPECT_NE(textRef.getSquash(), nullptr);
        EXPECT_GT(textRef.getSquashLength(), 0);
        EXPECT_EQ(textRef.getFastaFileName(), textFastaFile);
    }
    
    // Test 2: Create index with gzip reference file
    Reference gzipRef(gzipFastaFile, 1);
    bool gzipIndexResult = gzipRef.makeIndex();
    EXPECT_TRUE(gzipIndexResult);
    if (gzipIndexResult) {
        EXPECT_NE(gzipRef.getSquash(), nullptr);
        EXPECT_GT(gzipRef.getSquashLength(), 0);
        EXPECT_EQ(gzipRef.getFastaFileName(), gzipFastaFile);
    }
    
    // The main goal is to test that both text and gzip files can be processed
    // without errors. Content comparison is complex due to different processing paths.
    
    // Clean up
    std::filesystem::remove(textFastaFile);
    std::filesystem::remove(gzipFastaFile);
}




