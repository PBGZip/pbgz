/*
 * bam_record_testcase.cpp - Tests for the SAM line -> BAM record conversion
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

#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "bam_record.h"

/*
 * Two fields of a written record cannot be checked by reading the file back as SAM, because SAM
 * does not carry them: `bin` (index lookup) and the type letter of an integer optional field
 * (SAM prints every integer subtype as 'i'). BAM is what a round-tripped file is compared
 * against, so both are pinned here to the values htslib writes for the same SAM line.
 */
namespace {

typedef std::map<std::string, int32_t> RefIndex;

RefIndex makeRefIndex()
{
    RefIndex idx;
    idx["chr1"] = 0;
    idx["chr2"] = 1;
    return idx;
}

/* The record's bytes without the leading block_size. */
std::vector<uint8_t> build(const std::string& line, const RefIndex& idx)
{
    bamrec::BamRecordScratch scratch;
    EXPECT_EQ(1, bamrec::buildBamRecordFromSamLine((const uint8_t*)line.data(), line.size(),
                                                  idx, scratch));
    if (scratch.rec.size() < 4) {
        return std::vector<uint8_t>();
    }
    return std::vector<uint8_t>(scratch.rec.begin() + 4, scratch.rec.end());
}

uint16_t u16At(const std::vector<uint8_t>& rec, size_t off)
{
    uint16_t v = 0;
    std::memcpy(&v, rec.data() + off, 2);
    return v;
}

int32_t i32At(const std::vector<uint8_t>& rec, size_t off)
{
    int32_t v = 0;
    std::memcpy(&v, rec.data() + off, 4);
    return v;
}

/* bin sits after refID(4) + pos(4) + l_read_name(1) + mapq(1). */
uint16_t binOf(const std::vector<uint8_t>& rec)
{
    return u16At(rec, 10);
}

/* Optional fields start after read_name, CIGAR, SEQ (4-bit packed) and QUAL. */
size_t auxOffset(const std::vector<uint8_t>& rec)
{
    const uint32_t lReadName = rec[8];
    const uint32_t nCigar = u16At(rec, 12);
    const int32_t lSeq = i32At(rec, 16);
    return 32 + lReadName + nCigar * 4 + (size_t)((lSeq + 1) / 2) + (size_t)lSeq;
}

std::vector<uint8_t> auxBytes(const std::vector<uint8_t>& rec)
{
    const size_t off = auxOffset(rec);
    return std::vector<uint8_t>(rec.begin() + off, rec.end());
}

std::string samLine(const std::string& qname, uint16_t flag, const std::string& rname,
                    int64_t pos, const std::string& cigar, const std::string& seq,
                    const std::string& opt)
{
    const std::string qual(seq == "*" ? 1 : seq.size(), 'I');
    std::string s = qname + "\t" + std::to_string(flag) + "\t" + rname + "\t" +
                    std::to_string(pos) + "\t60\t" + cigar + "\t*\t0\t0\t" + seq + "\t" +
                    (seq == "*" ? "*" : qual);
    if (!opt.empty()) {
        s += "\t" + opt;
    }
    return s;
}

/* reg2bin of the span, as htslib computes it for a read that has a coordinate. */
TEST(BamRecordTest, MappedReadBinIsReg2BinOfTheSpan)
{
    const RefIndex idx = makeRefIndex();
    /* chr1:1 + 100M -> a 100 bp span wholly inside the first 16 kb bin. */
    EXPECT_EQ(4681, binOf(build(samLine("m", 0, "chr1", 1, "100M", std::string(100, 'A'), ""), idx)));
    /* 1M20000D1M spans 20001 bp, past the 16 kb bin: 585 from samtools for the same line. */
    EXPECT_EQ(585, binOf(build(samLine("d", 0, "chr1", 1, "1M20000D1M", "AC", ""), idx)));
    /* 1M2000000D1M spans 2 Mb, which reaches the top-level bin: 9. */
    EXPECT_EQ(9, binOf(build(samLine("b", 0, "chr1", 1, "1M2000000D1M", "AC", ""), idx)));
    /* A mapped read with no CIGAR is one reference base long. */
    EXPECT_EQ(4681, binOf(build(samLine("n", 0, "chr1", 1, "*", "AC", ""), idx)));
}

TEST(BamRecordTest, RecordWithoutCoordinateGetsTheUnmappedBin)
{
    const RefIndex idx = makeRefIndex();
    /* RNAME '*' and POS 0: unmapped (t2k's whole file is these). */
    EXPECT_EQ(bamrec::kUnmappedBin,
              binOf(build(samLine("u", 4, "*", 0, "*", "ACGT", ""), idx)));
    /* POS 0 with a reference name: htslib treats the record as unmapped ("mapped query cannot
       have zero coordinate"), and its bin goes with that. */
    EXPECT_EQ(bamrec::kUnmappedBin,
              binOf(build(samLine("z", 0, "chr1", 0, "100M", std::string(100, 'A'), ""), idx)));
    /* An unmapped record that still carries a coordinate keeps the coordinate's bin, as htslib's
       does - the flag is not what decides it. */
    EXPECT_EQ(4681, binOf(build(samLine("f", 4, "chr1", 1, "*", "ACGT", ""), idx)));
}

/*
 * The type letter of an integer optional field: htslib takes the narrowest width and, for a
 * non-negative value, the unsigned letter. Writing the signed letter first made every NM:C:1
 * come back as NM:c:1 - the same value, a different byte, and a different file than the one
 * samtools writes.
 */
TEST(BamRecordTest, IntegerOptionalFieldsUseHtslibTypeLetters)
{
    const RefIndex idx = makeRefIndex();
    const std::vector<uint8_t> aux = auxBytes(build(
        samLine("a", 0, "chr1", 1, "10M", std::string(10, 'A'),
                "NM:i:1\tAS:i:300\tXS:i:70000\tXD:i:-5\tXE:i:-40000"), idx));

    const uint8_t expected[] = {
        'N', 'M', 'C', 1,                                   /* 1      -> uint8  */
        'A', 'S', 'S', 0x2C, 0x01,                          /* 300    -> uint16 */
        'X', 'S', 'I', 0x70, 0x11, 0x01, 0x00,              /* 70000  -> uint32 */
        'X', 'D', 'c', 0xFB,                                /* -5     -> int8   */
        'X', 'E', 'i', 0xC0, 0x63, 0xFF, 0xFF,              /* -40000 -> int32  */
    };
    ASSERT_EQ(sizeof(expected), aux.size());
    EXPECT_EQ(0, std::memcmp(expected, aux.data(), sizeof(expected)));
}

/* The boundary values, where the width and therefore the letter changes. */
TEST(BamRecordTest, IntegerOptionalFieldWidthBoundaries)
{
    const RefIndex idx = makeRefIndex();
    struct Case {
        const char* value;
        uint8_t type;
        size_t width;
    };
    const Case cases[] = {
        {"0",          'C', 1},
        {"255",        'C', 1},
        {"256",        'S', 2},
        {"65535",      'S', 2},
        {"65536",      'I', 4},
        {"-1",         'c', 1},
        {"-128",       'c', 1},
        {"-129",       's', 2},
        {"-32768",     's', 2},
        {"-32769",     'i', 4},
    };
    for (const Case& c : cases) {
        const std::vector<uint8_t> aux = auxBytes(
            build(samLine("a", 0, "chr1", 1, "10M", std::string(10, 'A'),
                          std::string("ZZ:i:") + c.value), idx));
        ASSERT_EQ(3 + c.width, aux.size()) << "value " << c.value;
        EXPECT_EQ('Z', aux[0]);
        EXPECT_EQ('Z', aux[1]);
        EXPECT_EQ(c.type, aux[2]) << "value " << c.value;
    }
}

/* Non-integer optional fields carry their type through untouched. */
TEST(BamRecordTest, OtherOptionalFieldTypesArePreserved)
{
    const RefIndex idx = makeRefIndex();
    const std::vector<uint8_t> aux = auxBytes(build(
        samLine("a", 0, "chr1", 1, "10M", std::string(10, 'A'),
                "XA:A:a\tXZ:Z:hi\tXH:H:0A0B\tXF:f:1.5\tXB:B:i,1,2"), idx));

    std::vector<uint8_t> expected;
    const uint8_t head[] = {
        'X', 'A', 'A', 'a',
        'X', 'Z', 'Z', 'h', 'i', 0,
        'X', 'H', 'H', '0', 'A', '0', 'B', 0,
        'X', 'F', 'f',
    };
    expected.insert(expected.end(), head, head + sizeof(head));
    const float f = 1.5f;
    const uint8_t* fp = (const uint8_t*)&f;
    expected.insert(expected.end(), fp, fp + 4);
    const uint8_t tail[] = {'X', 'B', 'B', 'i', 2, 0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0};
    expected.insert(expected.end(), tail, tail + sizeof(tail));

    ASSERT_EQ(expected.size(), aux.size());
    EXPECT_EQ(0, std::memcmp(expected.data(), aux.data(), expected.size()));
}

}  // namespace
