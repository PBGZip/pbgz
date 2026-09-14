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

#include <gtest/gtest.h>

#include <cstdint>
#include <cstring>
#include <map>
#include <string>
#include <vector>

#include "sam_field_layout.h"
#include "sam_field_rules.h"

namespace {

/* A TLEN exception stream in the layout this build writes: (delta line index, zigzag32 value)
   varint pairs, the line indices cumulative and strictly increasing. */
std::vector<uint8_t> tlenVarintPayload(const std::vector<std::pair<uint32_t, int32_t>>& entries)
{
    std::vector<uint8_t> out;
    uint32_t prev = 0;
    for (const auto& entry : entries) {
        appendVarint(out, entry.first - prev);
        prev = entry.first;
        appendVarint(out, tlenZigzag32(entry.second));
    }
    return out;
}

/* The same entries in the layout up to 2026-08-19: (line, value) as fixed int32 pairs, the line an
   absolute content index and the value a plain signed integer. */
std::vector<uint8_t> tlenFixedPayload(const std::vector<std::pair<uint32_t, int32_t>>& entries)
{
    std::vector<uint8_t> out(entries.size() * 2 * sizeof(uint32_t));
    for (size_t i = 0; i < entries.size(); ++i) {
        const uint32_t line = entries[i].first;
        const int32_t value = entries[i].second;
        memcpy(out.data() + i * 8, &line, sizeof(uint32_t));
        memcpy(out.data() + i * 8 + 4, &value, sizeof(uint32_t));
    }
    return out;
}

std::vector<uint8_t> absolutePositions(const std::vector<uint32_t>& positions)
{
    std::vector<uint8_t> out(positions.size() * sizeof(uint32_t));
    if (!positions.empty()) {
        memcpy(out.data(), positions.data(), out.size());
    }
    return out;
}

/* A sub-stream meta as the decoder sees it, built the way the encoder builds it: a name, and the
   keys the layout in question carries. */
Json::Value streamMeta(const char* sname)
{
    Json::Value meta(Json::objectValue);
    meta["sname"] = sname;
    return meta;
}

}  // namespace

TEST(SamFieldLayoutTest, TlenVarintLayoutDecodes) {
    /* Deltas stay small, and a negative TLEN (the second read of a pair) survives zigzag. */
    const std::vector<std::pair<uint32_t, int32_t>> entries = {{3, 276}, {7, -276}, {8, 0}, {1000, 12345}};
    const std::vector<uint8_t> payload = tlenVarintPayload(entries);

    std::map<uint32_t, int32_t> decoded;
    ASSERT_TRUE(tlenDecodeVarints(payload.data(), (uint32_t)payload.size(), 4, decoded));
    ASSERT_EQ(decoded.size(), entries.size());
    for (const auto& entry : entries) {
        ASSERT_TRUE(decoded.count(entry.first)) << "line " << entry.first;
        EXPECT_EQ(decoded[entry.first], entry.second) << "line " << entry.first;
    }
}

TEST(SamFieldLayoutTest, TlenVarintLayoutRejectsWhatItCannotAccountFor) {
    const std::vector<std::pair<uint32_t, int32_t>> entries = {{3, 276}, {7, -276}};
    const std::vector<uint8_t> payload = tlenVarintPayload(entries);
    std::map<uint32_t, int32_t> decoded;

    /* The count the field declares is what the reader checks itself against: a stream carrying a
       different number of entries is not this layout, whatever the bytes look like. */
    EXPECT_FALSE(tlenDecodeVarints(payload.data(), (uint32_t)payload.size(), 3, decoded));
    EXPECT_FALSE(tlenDecodeVarints(payload.data(), (uint32_t)payload.size(), 1, decoded));
    /* Truncation is rejected, both mid-entry and mid-varint. */
    EXPECT_FALSE(tlenDecodeVarints(payload.data(), (uint32_t)payload.size() - 1, 2, decoded));
    EXPECT_FALSE(tlenDecodeVarints(payload.data(), 1, 2, decoded));
    /* A varint longer than five bytes cannot come from a 32-bit value. */
    const std::vector<uint8_t> overlong = {0x80, 0x80, 0x80, 0x80, 0x80, 0x00, 0x00};
    EXPECT_FALSE(tlenDecodeVarints(overlong.data(), (uint32_t)overlong.size(), 0, decoded));
}

/*
 * The layout up to 2026-08-19, which an archive with no "exc" marker holds: (line, TLEN) as fixed
 * int32 pairs, eight bytes per entry. Nothing about it is self-describing beyond its size, which is
 * why the reader checks that size against the count the field declares.
 */
TEST(SamFieldLayoutTest, TlenFixedLayoutDecodes) {
    const std::vector<std::pair<uint32_t, int32_t>> entries = {{0, 123}, {5, -276}, {900, 0}};
    const std::vector<uint8_t> payload = tlenFixedPayload(entries);

    std::map<uint32_t, int32_t> decoded;
    ASSERT_TRUE(tlenDecodeFixedPairs(payload.data(), (uint32_t)payload.size(), 3, decoded));
    ASSERT_EQ(decoded.size(), entries.size());
    for (const auto& entry : entries) {
        ASSERT_TRUE(decoded.count(entry.first)) << "line " << entry.first;
        EXPECT_EQ(decoded[entry.first], entry.second) << "line " << entry.first;
    }
}

TEST(SamFieldLayoutTest, TlenFixedLayoutRejectsSizesThatDoNotAddUp) {
    const std::vector<uint8_t> payload = tlenFixedPayload({{0, 123}, {5, -276}});
    std::map<uint32_t, int32_t> decoded;

    /* Eight bytes per entry is the one thing this layout must satisfy. */
    EXPECT_FALSE(tlenDecodeFixedPairs(payload.data(), (uint32_t)payload.size(), 1, decoded));
    EXPECT_FALSE(tlenDecodeFixedPairs(payload.data(), (uint32_t)payload.size(), 3, decoded));
    EXPECT_FALSE(tlenDecodeFixedPairs(payload.data(), (uint32_t)payload.size() - 4, 2, decoded));
    EXPECT_FALSE(tlenDecodeFixedPairs(payload.data(), (uint32_t)payload.size(), 0, decoded));
}

/*
 * Why an archive without the marker can fall back on the size test: the fixed pairs do not survive
 * the varint reader. Eight bytes per entry make it walk off both ends of the entry count - it finds
 * roughly four times as many entries as were written, and consumes the buffer only by accident - so
 * the count the field declares is never met and the reader declines.
 */
TEST(SamFieldLayoutTest, TlenFixedBytesAreNotVarintPairs) {
    const std::vector<std::pair<uint32_t, int32_t>> entries = {{0, 123}, {5, -276}, {900, 4096}};
    const std::vector<uint8_t> payload = tlenFixedPayload(entries);

    std::map<uint32_t, int32_t> decoded;
    EXPECT_FALSE(tlenDecodeVarints(payload.data(), (uint32_t)payload.size(), 3, decoded));
}

/*
 * The first layout's SEQ exception list: a sub-stream named "npos" with no "ch" at all, holding
 * every position that layout knew of - 'N' and 'n' alike - with its count in the field-level
 * "ncount" and 'N' written back at each position. This is the one older SEQ layout still read.
 */
TEST(SamFieldLayoutTest, SeqLegacyNposIsReadAsUpperN) {
    Json::Value fieldMeta(Json::objectValue);
    fieldMeta["ncount"] = (Json::UInt)3;
    const Json::Value stream = streamMeta("npos");

    const SeqExcLayout layout = seqExcLayoutOf(fieldMeta, stream);
    EXPECT_EQ(layout.form, SeqExcForm::Absolute);
    EXPECT_EQ(layout.ch, 'N');
    EXPECT_EQ(layout.count, 3u);
    EXPECT_EQ(layout.runs, 0u);
}

TEST(SamFieldLayoutTest, SeqRoutingReadsTheNameAndTheStreamOwnKeys) {
    Json::Value fieldMeta(Json::objectValue);
    fieldMeta["ncount"] = (Json::UInt)99;   /* must not be consulted when the stream carries a count */

    /* The delta form: the character and the count travel in the stream's own meta. */
    Json::Value delta = streamMeta("nposd");
    delta["ch"] = (Json::UInt)'R';
    delta["count"] = (Json::UInt)2;
    const SeqExcLayout deltaLayout = seqExcLayoutOf(fieldMeta, delta);
    EXPECT_EQ(deltaLayout.form, SeqExcForm::Delta);
    EXPECT_EQ(deltaLayout.ch, 'R');
    EXPECT_EQ(deltaLayout.count, 2u);

    /* The run form adds how many runs the payload describes. */
    Json::Value run = streamMeta("nposr");
    run["ch"] = (Json::UInt)'n';
    run["count"] = (Json::UInt)6;
    run["runs"] = (Json::UInt)2;
    const SeqExcLayout runLayout = seqExcLayoutOf(fieldMeta, run);
    EXPECT_EQ(runLayout.form, SeqExcForm::Run);
    EXPECT_EQ(runLayout.ch, 'n');
    EXPECT_EQ(runLayout.count, 6u);
    EXPECT_EQ(runLayout.runs, 2u);

    /* The absolute form as it is written now: the same shape as the legacy one, plus "ch". */
    Json::Value abs = streamMeta("npos");
    abs["ch"] = (Json::UInt)'Y';
    abs["count"] = (Json::UInt)1;
    const SeqExcLayout absLayout = seqExcLayoutOf(fieldMeta, abs);
    EXPECT_EQ(absLayout.form, SeqExcForm::Absolute);
    EXPECT_EQ(absLayout.ch, 'Y');
    EXPECT_EQ(absLayout.count, 1u);
}

/*
 * Only the three forms are exception streams. Any other name - a sub-stream of another kind, or the
 * middle layout's "nposl"/"nposx", which this build reads no longer - ends the run of exception
 * streams and is never guessed at.
 */
TEST(SamFieldLayoutTest, SeqRoutingStopsAtAnyOtherName) {
    EXPECT_EQ(seqExcFormOf("nposd"), SeqExcForm::Delta);
    EXPECT_EQ(seqExcFormOf("npos"), SeqExcForm::Absolute);
    EXPECT_EQ(seqExcFormOf("nposr"), SeqExcForm::Run);
    EXPECT_EQ(seqExcFormOf("baselen"), SeqExcForm::None);
    EXPECT_EQ(seqExcFormOf("nposl"), SeqExcForm::None);
    EXPECT_EQ(seqExcFormOf("nposx"), SeqExcForm::None);
    EXPECT_EQ(seqExcFormOf(""), SeqExcForm::None);
}

TEST(SamFieldLayoutTest, SeqAbsolutePayloadIsOneOffsetPerPosition) {
    Json::Value fieldMeta(Json::objectValue);
    fieldMeta["ncount"] = (Json::UInt)3;
    const SeqExcLayout layout = seqExcLayoutOf(fieldMeta, streamMeta("npos"));
    const std::vector<uint8_t> payload = absolutePositions({5, 40, 4096});

    std::vector<uint32_t> pos;
    ASSERT_TRUE(seqExcDecodePositions(layout, payload.data(), (uint32_t)payload.size(), pos));
    EXPECT_EQ(pos, std::vector<uint32_t>({5, 40, 4096}));

    /* The size has to be the count times four: anything else is not this layout. */
    EXPECT_FALSE(seqExcDecodePositions(layout, payload.data(), (uint32_t)payload.size() - 1, pos));
    EXPECT_FALSE(seqExcDecodePositions(layout, payload.data(), (uint32_t)payload.size() + 4, pos));
}

TEST(SamFieldLayoutTest, SeqDeltaPayloadIsCumulative) {
    Json::Value fieldMeta(Json::objectValue);
    Json::Value stream = streamMeta("nposd");
    stream["ch"] = (Json::UInt)'N';
    stream["count"] = (Json::UInt)4;
    const SeqExcLayout layout = seqExcLayoutOf(fieldMeta, stream);

    std::vector<uint8_t> payload;
    for (uint32_t delta : {3u, 4u, 100u, 1u}) {
        appendVarint(payload, delta);
    }

    std::vector<uint32_t> pos;
    ASSERT_TRUE(seqExcDecodePositions(layout, payload.data(), (uint32_t)payload.size(), pos));
    EXPECT_EQ(pos, std::vector<uint32_t>({3, 7, 107, 108}));
    /* One entry short is a truncated stream, not a shorter list. */
    EXPECT_FALSE(seqExcDecodePositions(layout, payload.data(), (uint32_t)payload.size() - 1, pos));
}

TEST(SamFieldLayoutTest, SeqRunPayloadIsGapsThenLengths) {
    Json::Value fieldMeta(Json::objectValue);
    Json::Value stream = streamMeta("nposr");
    stream["ch"] = (Json::UInt)'N';
    stream["count"] = (Json::UInt)5;
    stream["runs"] = (Json::UInt)2;
    const SeqExcLayout layout = seqExcLayoutOf(fieldMeta, stream);

    /* Two runs: the first starts 3 in and covers 3 positions (3,4,5); the second starts 2 after
       that run's end (6) and covers 2 (8,9). A gap is measured from the end of the previous run,
       which is what the encoder writes (see SeqExceptionClass::closeRun). */
    std::vector<uint8_t> payload;
    appendVarint(payload, 3);
    appendVarint(payload, 2);
    appendVarint(payload, 3);
    appendVarint(payload, 2);

    std::vector<uint32_t> pos;
    ASSERT_TRUE(seqExcDecodePositions(layout, payload.data(), (uint32_t)payload.size(), pos));
    EXPECT_EQ(pos, std::vector<uint32_t>({3, 4, 5, 8, 9}));

    /* The runs have to describe exactly the count the field declares. */
    stream["count"] = (Json::UInt)4;
    const SeqExcLayout shortLayout = seqExcLayoutOf(fieldMeta, stream);
    EXPECT_FALSE(seqExcDecodePositions(shortLayout, payload.data(), (uint32_t)payload.size(), pos));
}

/*
 * The template length both actuators encode TLEN against: the two aligner conventions differ by
 * exactly 1, and the sign follows the side the record is on.
 */
TEST(SamFieldRuleTest, TemplateLengthUnderBothConventions) {
    /* Same pair, 76-base reads, left end at 100 and right end at 300. */
    EXPECT_EQ(samTemplateLen(100, 300, 76, 76, false), 276);
    EXPECT_EQ(samTemplateLen(100, 300, 76, 76, true), 275);
    /* The right end reports the same template with the opposite sign. */
    EXPECT_EQ(samTemplateLen(300, 100, 76, 76, false), -276);
    EXPECT_EQ(samTemplateLen(300, 100, 76, 76, true), -275);
    /* A mate whose span is not in this block contributes 0, which only costs ratio. */
    EXPECT_EQ(samTemplateLen(100, 300, 76, 0, false), 200);
    /* Both ends at the same position: the mate's side decides, as in the positive branch. */
    EXPECT_EQ(samTemplateLen(300, 300, 76, 50, false), 50);
}
