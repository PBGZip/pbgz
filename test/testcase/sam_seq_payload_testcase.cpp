/*
 * sam_field_layout_testcase.cpp - the archive's SAM field payload layouts, current and older.
 * Copyright (C) 2025 PBGZip
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/*
 * The readers in sam_field_layout.h exist for archives older than this build, and an archive from a
 * released revision is not something a test can be written against - the point of the module is
 * that a layout is a pure function of a payload, so each one is pinned here by its bytes: the
 * layouts this build writes, the layouts older revisions wrote, and the cases that must be
 * rejected rather than guessed at.
 *
 * The two that matter most for compatibility are the TLEN exception stream of the revisions up to
 * 2026-08-19 (fixed int32 pairs, no marker - the layout an archive with no "exc" key holds) and the
 * first layout's lone SEQ "npos" (absolute offsets, no "ch", counted by the field-level "ncount",
 * 'N' and 'n' alike). Both are exercised below, together with the bytes that tell them apart from
 * the layouts written now.
 */

/*
 * sam_seq_payload_testcase.cpp - the shared SEQ payload machinery: the CIGAR parse both spans come
 * out of, and the walk the writer and the reader put those spans through.
 *
 * These live apart from the field-layout cases because they are the module's own contract rather
 * than an archive layout: the spans and the walk decide what the payloads mean, and both sides of
 * the codec hang off them.
 */

#include <gtest/gtest.h>

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

#include "sam_info.h"
#include "sam_seq_payload.h"

/*
 * The CIGAR walk both SEQ paths run (see seqWalkCigar): what it reports for a read the ops describe,
 * what it reports when they run past it, and that H/P contribute nothing while D/N move only the
 * reference. The walk is the shared piece - the writer's fallback and the reader's reconstruction
 * both hang off its answer - so it is pinned here rather than through either of them.
 */
TEST(SamSeqPayloadTest, CigarWalkReportsSegmentsAndReferenceSpan) {
    const std::vector<CigarOp> ops = {{'M', 4}, {'I', 2}, {'D', 3}, {'M', 4}, {'H', 5}, {'S', 2}};
    std::string seen;
    const SeqWalkEnd end = seqWalkCigar(ops, 12, 100,
        [&seen](uint32_t readPos, int64_t refPos, uint32_t len) {
            seen += "M" + std::to_string(readPos) + "@" + std::to_string(refPos) + "x" +
                    std::to_string(len) + ";";
            return SeqWalkAction::Continue;
        },
        [&seen](uint32_t readPos, uint32_t len) {
            seen += "I" + std::to_string(readPos) + "x" + std::to_string(len) + ";";
            return SeqWalkAction::Continue;
        },
        [&seen](uint32_t len) {
            seen += "D" + std::to_string(len) + ";";
            return SeqWalkAction::Continue;
        });
    EXPECT_TRUE(end.complete);
    EXPECT_EQ(end.readPos, 12u);
    EXPECT_EQ(end.refPos, 111);   /* 100, +4 for the first M, +3 for the D, +4 for the second */
    EXPECT_EQ(seen, std::string("M0@100x4;I4x2;D3;M6@107x4;I10x2;"));
}

TEST(SamSeqPayloadTest, CigarWalkEndsWhenOpsRunPastTheRead) {
    const std::vector<CigarOp> ops = {{'M', 4}, {'M', 6}};
    uint32_t handed = 0;
    const SeqWalkEnd end = seqWalkCigar(ops, 8, 0,
        [&handed](uint32_t, int64_t, uint32_t) { ++handed; return SeqWalkAction::Continue; },
        [](uint32_t, uint32_t) { return SeqWalkAction::Continue; },
        [](uint32_t) { return SeqWalkAction::Continue; });
    EXPECT_FALSE(end.complete);   /* 4 + 6 > 8, so the ops do not describe the read */
    EXPECT_EQ(end.readPos, 4u);
    EXPECT_EQ(handed, 1u);        /* the op that overruns is never handed over */
}

TEST(SamSeqPayloadTest, CigarWalkWithoutOpsIsIncompleteForANonEmptyRead) {
    const std::vector<CigarOp> ops;
    uint32_t handed = 0;
    const SeqWalkEnd end = seqWalkCigar(ops, 3, 7,
        [&handed](uint32_t, int64_t, uint32_t) { ++handed; return SeqWalkAction::Continue; },
        [](uint32_t, uint32_t) { return SeqWalkAction::Continue; },
        [](uint32_t) { return SeqWalkAction::Continue; });
    EXPECT_FALSE(end.complete);
    EXPECT_EQ(handed, 0u);
    EXPECT_EQ(end.refPos, 7);
}

/*
 * The one CIGAR parse (see sam_seq_payload.h): the two spans it answers and the operation list. The
 * lower-case forms are pinned as they are - the read span counts a lower-case 'x', the reference
 * span has never claimed to - because both spans feed what the archive records, so the difference is
 * a decision rather than an accident to be tidied away.
 */
TEST(SamSeqPayloadTest, CigarParseAnswersBothSpansAndTheOps) {
    const char* cigar = "6S30M1I114S";
    const uint32_t len = (uint32_t)std::strlen(cigar);
    std::vector<CigarOp> ops;
    EXPECT_EQ(parseCigarOps((const uint8_t*)cigar, len, ops), 30u);   /* only the M */
    EXPECT_EQ(cigarSeqConsumed((const uint8_t*)cigar, len), 151u);    /* 6 + 30 + 1 + 114 */
    ASSERT_EQ(ops.size(), 4u);
    EXPECT_EQ(ops[0].op, 'S');
    EXPECT_EQ(ops[0].len, 6u);
    EXPECT_EQ(ops[3].op, 'S');
    EXPECT_EQ(ops[3].len, 114u);
}

TEST(SamSeqPayloadTest, CigarParseKeepsEachSpansLetterSet) {
    /* H and P count towards neither span; D/N towards the reference span only; I/S towards the read
       span only; the operation list holds the upper-case forms alone. */
    const char* cigar = "10H5D3I2P7N4S1X2x";
    const uint32_t len = (uint32_t)std::strlen(cigar);
    std::vector<CigarOp> ops;
    parseCigarOps((const uint8_t*)cigar, len, ops);
    EXPECT_EQ(cigarRefConsumed((const uint8_t*)cigar, len), 5u + 7u + 1u);
    EXPECT_EQ(cigarSeqConsumed((const uint8_t*)cigar, len), 3u + 4u + 1u + 2u);
    /* The list holds the upper-case operations, H and P included: H, D, I, P, N, S and X. */
    ASSERT_EQ(ops.size(), 7u);
    EXPECT_EQ(ops[0].op, 'H');
    EXPECT_EQ(ops[6].op, 'X');
}

/*
 * referenceLinesUpWithHeader is the gate for the reference-coded SEQ layout. It has to accept a
 * FASTA of the same assembly even when its contig names follow the other convention (Ensembl
 * "1".."29" over a RefSeq "NC_037328.1"... FASTA), and it has to reject a reference of a different
 * assembly. The alias case is decided by position because that is how SEQ addresses the reference:
 * every @SQ contig must start at the offset the reference gives it.
 */
TEST(SamSeqPayloadTest, ReferenceLinesUpWithHeaderAcceptsAliasNamesByPosition)
{
    const std::string fasta = "test_alias_ref.fa";
    {
        std::ofstream out(fasta);
        ASSERT_TRUE(out.is_open());
        out << ">NC_000001.1\n" << std::string(40, 'A') << '\n';
        out << ">NC_000002.1\n" << std::string(24, 'C') << '\n';
    }
    Reference reference(fasta, 1);
    ASSERT_TRUE(reference.makeIndex());

    SamInfo& info = SamInfo::getInstance();
    info.clearChromosomeInfo();
    info.resetChrIdCounter();
    info.addChromosomeInfo("1", 40);
    info.addChromosomeInfo("2", 24);
    EXPECT_TRUE(referenceLinesUpWithHeader(&reference));

    info.clearChromosomeInfo();
    info.resetChrIdCounter();
    std::remove(fasta.c_str());
}

TEST(SamSeqPayloadTest, ReferenceLinesUpWithHeaderAcceptsCanonicalNamesInAnyOrder)
{
    const std::string fasta = "test_canonical_ref.fa";
    {
        std::ofstream out(fasta);
        ASSERT_TRUE(out.is_open());
        /* "chr1"/"chr2" against "1"/"2", in the opposite order: only the canonical name matches. */
        out << ">2\n" << std::string(24, 'C') << '\n';
        out << ">1\n" << std::string(40, 'A') << '\n';
    }
    Reference reference(fasta, 1);
    ASSERT_TRUE(reference.makeIndex());

    SamInfo& info = SamInfo::getInstance();
    info.clearChromosomeInfo();
    info.resetChrIdCounter();
    info.addChromosomeInfo("chr1", 40);
    info.addChromosomeInfo("chr2", 24);
    EXPECT_TRUE(referenceLinesUpWithHeader(&reference));

    info.clearChromosomeInfo();
    info.resetChrIdCounter();
    std::remove(fasta.c_str());
}

TEST(SamSeqPayloadTest, ReferenceLinesUpWithHeaderRejectsAnotherAssembly)
{
    const std::string fasta = "test_other_assembly_ref.fa";
    {
        std::ofstream out(fasta);
        ASSERT_TRUE(out.is_open());
        /* Same names, different lengths and layout: another assembly. */
        out << ">NC_000001.1\n" << std::string(48, 'A') << '\n';
        out << ">NC_000002.1\n" << std::string(24, 'C') << '\n';
    }
    Reference reference(fasta, 1);
    ASSERT_TRUE(reference.makeIndex());

    SamInfo& info = SamInfo::getInstance();
    info.clearChromosomeInfo();
    info.resetChrIdCounter();
    info.addChromosomeInfo("1", 40);
    info.addChromosomeInfo("2", 24);
    EXPECT_FALSE(referenceLinesUpWithHeader(&reference));

    info.clearChromosomeInfo();
    info.resetChrIdCounter();
    std::remove(fasta.c_str());
}

/*
 * The split SEQ layout separates aligned bytes (M/=/X) from inserted/soft-clipped ones, so the
 * decoder can use a different model for each. Partitioning must follow the CIGAR and keep read
 * order; a record that was not reference-coded belongs entirely to the inserted/raw stream.
 */
TEST(SamSeqPayloadTest, PartitionSeqPayloadSplitsAlignedAndInsertedBytes)
{
    const std::vector<CigarOp> ops = {{'M', 4}, {'I', 2}, {'M', 3}, {'S', 2}, {'D', 1}, {'M', 1}};
    const uint8_t coded[12] = {0, 1, 2, 3, 10, 11, 4, 5, 6, 12, 13, 7};
    std::vector<uint8_t> mat;
    std::vector<uint8_t> ins;
    ASSERT_TRUE(partitionSeqPayload(ops, coded, 12, true, mat, ins));
    EXPECT_EQ(mat, (std::vector<uint8_t>{0, 1, 2, 3, 4, 5, 6, 7}));
    EXPECT_EQ(ins, (std::vector<uint8_t>{10, 11, 12, 13}));
}

TEST(SamSeqPayloadTest, PartitionSeqPayloadSendsRawRecordsToInsertedStream)
{
    const uint8_t raw[3] = {'A', 'C', 'G'};
    std::vector<uint8_t> mat;
    std::vector<uint8_t> ins;
    ASSERT_TRUE(partitionSeqPayload({}, raw, 3, false, mat, ins));
    EXPECT_TRUE(mat.empty());
    EXPECT_EQ(ins, (std::vector<uint8_t>{'A', 'C', 'G'}));
}

TEST(SamSeqPayloadTest, PartitionSeqPayloadFollowsTheReferenceWalk)
{
    /* D/N consume nothing, I/S go to the inserted stream, and H/P are ignored. */
    const std::vector<CigarOp> ops = {{'H', 2}, {'S', 1}, {'M', 2}, {'N', 3}, {'=', 1}, {'I', 1}};
    const uint8_t coded[5] = {9, 0, 1, 2, 3};
    std::vector<uint8_t> mat;
    std::vector<uint8_t> ins;
    ASSERT_TRUE(partitionSeqPayload(ops, coded, 5, true, mat, ins));
    EXPECT_EQ(mat, (std::vector<uint8_t>{0, 1, 2}));
    EXPECT_EQ(ins, (std::vector<uint8_t>{9, 3}));
}
