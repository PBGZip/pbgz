/*
 * coder_fcv2.cpp - Implementation of the quality-value context-mixing coder.
 *
 * A quality symbol is split by a Huffman tree into binary decisions; several
 * models of differing granularity predict at each decision and their log-odds
 * are mixed with gradient-descent-updated weights before range coding.
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

#include "coder_fcv2.h"

#include <string.h>
#include <math.h>

#include <limits>

#include "coder.h"
#include "fc/rangecoder.h"

namespace {

/* Upper bound on the number of Huffman tree nodes; also limits the alphabet size. */
const int TREE_CAP = FCV2_TREE_CAP;

/* Weight-update step size as a shift, equivalent to a learning rate of 1/4096. */
const int WEIGHT_LR_SHIFT = 12;

/* The mix is fixed at FCV2_MODEL_COUNT models; the arrays below are sized by it. */

/* Legal range of the context parameter tiers; the decoder normalizes by the same
   rule so a corrupted stream cannot allocate an oversized model array. Every
   maximum is written into the stream header as one byte. */
const int CFG_CYCLE_MAX_MIN = 32;
const int CFG_CYCLE_MAX_MAX = 128;
const int CFG_CYCLE_BUCKET_MAX = 32;
const int CFG_PREV_SHIFT_MAX = 3;

static_assert(CFG_CYCLE_MAX_MAX <= 255 && CFG_CYCLE_BUCKET_MAX <= 255 &&
              CFG_PREV_SHIFT_MAX <= 255 && FCV2_MODEL_COUNT <= 255 &&
              TREE_CAP <= 255,
              "the tier fields and the alphabet size are written into the stream header "
              "one byte each, so a maximum above 255 would be truncated on write and "
              "break decoding");

/* Normalizes cfg to the tiers actually used; both sides must call this or the
   context indices misalign. */
Fcv2Cfg normalizeCfg(Fcv2Cfg cfg)
{
    if (cfg.cycleMax < CFG_CYCLE_MAX_MIN) cfg.cycleMax = CFG_CYCLE_MAX_MIN;
    if (cfg.cycleMax > CFG_CYCLE_MAX_MAX) cfg.cycleMax = CFG_CYCLE_MAX_MAX;
    if (cfg.cycleBucket < 1) cfg.cycleBucket = 1;
    if (cfg.cycleBucket > CFG_CYCLE_BUCKET_MAX) cfg.cycleBucket = CFG_CYCLE_BUCKET_MAX;
    if (cfg.prevShift < 0) cfg.prevShift = 0;
    if (cfg.prevShift > CFG_PREV_SHIFT_MAX) cfg.prevShift = CFG_PREV_SHIFT_MAX;
    return cfg;
}

/* Number of distinct base-context states: ACGTN + unknown (seq empty or a
   non-base character). */
const int BASE_STATES = 6;

/* Raw byte -> base symbol (ACGTN -> 0..4), everything else maps to 5
   (unknown). Filled in initTables. */
uint8_t g_baseSym[256];

/* Model snapshot header: lets a snapshot be validated against the current array
   layout before it is read. */
const uint8_t MODEL_MAGIC[] = { 'f', 'c', 'v', '2', 'p', 'r', 'i', 'o', 'r' };
const uint16_t MODEL_FORMAT_VERSION = 9;  /* v9: m5/m6 removed, model count fixed at FCV2_MODEL_COUNT */

inline void appendU16(std::vector<uint8_t>& out, uint16_t value)
{
    out.push_back((uint8_t)(value & 0xFF));
    out.push_back((uint8_t)(value >> 8));
}

inline void appendU32(std::vector<uint8_t>& out, uint32_t value)
{
    for (int shift = 0; shift < 32; shift += 8) {
        out.push_back((uint8_t)(value >> shift));
    }
}

inline void appendI32(std::vector<uint8_t>& out, int32_t value)
{
    appendU32(out, (uint32_t)value);
}

/* Byte reader for the snapshot; advances only when enough input remains, so a
   truncated blob cannot read out of bounds. All numbers are little-endian. */
class ModelBlobReader {
public:
    explicit ModelBlobReader(const std::vector<uint8_t>& blob) : data(blob), pos(0) {}

    bool readU8(uint8_t& value)
    {
        if (pos == data.size()) return false;
        value = data[pos++];
        return true;
    }

    bool readU16(uint16_t& value)
    {
        uint8_t lo = 0, hi = 0;
        if (!readU8(lo) || !readU8(hi)) return false;
        value = (uint16_t)((uint16_t)lo | ((uint16_t)hi << 8));
        return true;
    }

    bool readU32(uint32_t& value)
    {
        value = 0;
        for (int shift = 0; shift < 32; shift += 8) {
            uint8_t byte = 0;
            if (!readU8(byte)) return false;
            value |= (uint32_t)byte << shift;
        }
        return true;
    }

    bool atEnd() const { return pos == data.size(); }

private:
    const std::vector<uint8_t>& data;
    size_t pos;
};

/* A probability and a count packed into one uint16_t: low 12 bits probability
   (0..4095), high 4 bits update count (0..15). */
typedef uint16_t Counter;

inline int counterProb(Counter c)  { return c & 0x0FFF; }
inline int counterCount(Counter c) { return c >> 12; }

const Counter COUNTER_INIT = 2048;   /* probability 0.5, count 0 */

/* Lookup tables converting between log-odds and probability, 12-bit fixed
   point with a scale factor of 256. */
short  g_stretch[4096];
uint16_t g_squash[4096];

/* Adaptive step-size table: large step when the update count is small, shrinking
   as it grows. */
uint8_t g_adaptShift[16];

bool g_tablesReady = false;

void initTables()
{
    if (g_tablesReady) {
        return;
    }
    for (int i = 0; i < 4096; i++) {
        double v = 4096.0 / (1.0 + exp(-(i - 2048) / 256.0));
        if (v < 1)    v = 1;
        if (v > 4095) v = 4095;
        g_squash[i] = (uint16_t)v;
    }
    for (int p = 0; p < 4096; p++) {
        int pp = p < 1 ? 1 : (p > 4095 ? 4095 : p);
        double v = 256.0 * log((double)pp / (4096.0 - pp));
        if (v < -2047) v = -2047;
        if (v >  2047) v =  2047;
        g_stretch[p] = (short)v;
    }
    for (int n = 0; n < 16; n++) {
        int s = 1 + n / 3;
        g_adaptShift[n] = (uint8_t)(s > 6 ? 6 : s);
    }
    for (int i = 0; i < 256; i++) g_baseSym[i] = BASE_STATES - 1;
    g_baseSym[(int)'A'] = g_baseSym[(int)'a'] = 0;
    g_baseSym[(int)'C'] = g_baseSym[(int)'c'] = 1;
    g_baseSym[(int)'G'] = g_baseSym[(int)'g'] = 2;
    g_baseSym[(int)'T'] = g_baseSym[(int)'t'] = 3;
    g_baseSym[(int)'N'] = g_baseSym[(int)'n'] = 4;
    g_tablesReady = true;
}

inline int squash(int x)
{
    if (x <= -2048) return 1;
    if (x >=  2047) return 4095;
    return g_squash[x + 2048];
}

inline void counterUpdate(Counter& c, int isZero)
{
    int target = isZero ? 4095 : 0;
    int p = counterProb(c);
    int n = counterCount(c);
    p += ((target - p) * 2) >> g_adaptShift[n];
    if (p < 1)    p = 1;
    if (p > 4095) p = 4095;
    if (n < 15)   n++;
    c = (Counter)(p | (n << 12));
}

/* Huffman tree: leaves are quality symbols; frequency-built so common symbols
   get shorter paths. Affects speed only, not the compression ratio. */
struct HuffTree {
    int alphaSize;
    int root;
    int child[TREE_CAP * 2][2];
    uint8_t pathLen[TREE_CAP];
    uint8_t pathNode[TREE_CAP][24];
    uint8_t pathBit[TREE_CAP][24];

    void build(const uint32_t* freq, int alpha)
    {
        alphaSize = alpha;
        long long weight[TREE_CAP * 2];
        int node[TREE_CAP * 2];
        int live = 0;
        for (int i = 0; i < alpha; i++) {
            /* Add one to the frequency so that even symbols never observed get
               a path, avoiding out-of-bounds access during coding. */
            weight[live] = (long long)freq[i] + 1;
            node[live] = i;
            live++;
        }
        int nextNode = alpha;
        while (live > 1) {
            int a = 0;
            for (int i = 1; i < live; i++) {
                if (weight[i] < weight[a]) a = i;
            }
            int b = -1;
            for (int i = 0; i < live; i++) {
                if (i != a && (b < 0 || weight[i] < weight[b])) b = i;
            }
            int parent = nextNode++;
            child[parent - alpha][0] = node[a];
            child[parent - alpha][1] = node[b];
            long long merged = weight[a] + weight[b];
            int keep = a < b ? a : b;
            int drop = a < b ? b : a;
            weight[keep] = merged;
            node[keep] = parent;
            for (int i = drop; i < live - 1; i++) {
                weight[i] = weight[i + 1];
                node[i] = node[i + 1];
            }
            live--;
        }
        root = node[0];
        memset(pathLen, 0, sizeof(pathLen));
        uint8_t stackNode[32];
        uint8_t stackBit[32];
        walk(root, 0, stackNode, stackBit);
    }

private:
    void walk(int n, int depth, uint8_t* sn, uint8_t* sb)
    {
        if (n < alphaSize) {
            pathLen[n] = (uint8_t)depth;
            for (int d = 0; d < depth; d++) {
                pathNode[n][d] = sn[d];
                pathBit[n][d] = sb[d];
            }
            return;
        }
        for (int b = 0; b < 2; b++) {
            sn[depth] = (uint8_t)(n - alphaSize);
            sb[depth] = (uint8_t)b;
            walk(child[n - alphaSize][b], depth + 1, sn, sb);
        }
    }
};

/*
 * Five context models of increasing granularity:
 *   m0 tree node; m1 position + node; m2 previous symbol + position + node;
 *   m3 previous two symbols + position bin + strand + node; m4 base + position + node.
 * Array sizes are determined by Fcv2Cfg; see init().
 */
struct ContextModel {
    int alphaSize;
    int prevStates;          /* number of symbol values plus one; the extra slot represents "no predecessor" */
    int prevQStates;         /* number of quantized predecessor values for m3 (>> prevShift) */
    int cycleMax;
    int cycleBucket;
    int prevShift;
    std::vector<Counter> m0, m1, m2, m3, m4;
    /* The five tables as base pointers, set once in init() after sizing. */
    Counter* table[FCV2_MODEL_COUNT] = { nullptr, nullptr, nullptr, nullptr, nullptr };

    /* Points `table` at the current buffers; called after init() and after a move. */
    void rebindTables()
    {
        table[0] = m0.data();
        table[1] = m1.data();
        table[2] = m2.data();
        table[3] = m3.data();
        table[4] = m4.data();
    }

    ContextModel() = default;
    ContextModel(const ContextModel&) = delete;
    ContextModel& operator=(const ContextModel&) = delete;
    ContextModel(ContextModel&&) = default;
    ContextModel& operator=(ContextModel&&) = default;
    std::vector<int> weight; /* mixing weights, organized by position bin x tree node x model */
    /* cycle index -> position bin, precomputed by init(); see there. */
    std::vector<int> cycBucketTab;

    void init(int alpha, const Fcv2Cfg& cfg)
    {
        alphaSize = alpha;
        prevStates = alpha + 1;
        cycleMax = cfg.cycleMax;
        cycleBucket = cfg.cycleBucket;
        prevShift = cfg.prevShift;
        /* Quantizing the previous quality value by >>prevShift keeps the order-2
           context from being too sparse. */
        prevQStates = (alpha >> prevShift) + 2;
        m0.assign((size_t)TREE_CAP, COUNTER_INIT);
        m1.assign((size_t)cycleMax * TREE_CAP, COUNTER_INIT);
        m2.assign((size_t)prevStates * cycleMax * TREE_CAP, COUNTER_INIT);
        m3.assign((size_t)prevQStates * prevQStates * cycleBucket * 2 * TREE_CAP, COUNTER_INIT);
        m4.assign((size_t)BASE_STATES * cycleMax * TREE_CAP, COUNTER_INIT);
        /* The mix loop's fast path: see the note on `table`. */
        rebindTables();
        /* Initial weights 1<<14 (fixed-point 0.25) start close to a simple average. */
        weight.assign((size_t)cycleBucket * TREE_CAP * FCV2_MODEL_COUNT, 1 << 14);
        /* Precomputed cycle-index -> position-bin map, shared by m3 and the
           weight row. */
        cycBucketTab.assign((size_t)cycleMax, 0);
        for (int c = 0; c < cycleMax; c++) {
            int b = c * cycleBucket / cycleMax;
            cycBucketTab[c] = (b >= cycleBucket) ? cycleBucket - 1 : b;
        }
    }

    inline void slotBases(int prev, int prev2, int cyc, int rev, int baseSym,
                          size_t* base, int& wBucket) const
    {
        int bucket = cycBucketTab[cyc];
        wBucket = bucket;   /* handed back to the caller: the weight row shares it */
        int qp = (prev >> prevShift);      /* quantized previous quality value */
        int qp2 = (prev2 >> prevShift);
        if (qp >= prevQStates) qp = prevQStates - 1;
        if (qp2 >= prevQStates) qp2 = prevQStates - 1;
        base[0] = 0;
        base[1] = (size_t)cyc * TREE_CAP;
        base[2] = ((size_t)prev * cycleMax + cyc) * TREE_CAP;
        base[3] = ((((size_t)qp2 * prevQStates + qp) * cycleBucket + bucket) * 2 + rev) * TREE_CAP;
        base[4] = ((size_t)baseSym * cycleMax + cyc) * TREE_CAP;
    }

    inline int* weightAt(int bucket, int node)
    {
        return &weight[((size_t)bucket * TREE_CAP + node) * FCV2_MODEL_COUNT];
    }

    inline Counter* counterAt(int model, size_t base, int node)
    {
        return table[model] + base + node;
    }
};

} /* namespace */

/* The implementation body; outside the anonymous namespace because the header
   forward-declares it. */
class fcv2_impl {
public:
    coder_io* io;
    ContextModel cm;
    HuffTree tree;
    FcRangeCoder rc;

    /*
     * Probability model for the strand. The strand is written into the stream
     * with each record so the coder is self-contained; the cost (about one bit
     * per record) is less than the gain it enables.
     */
    Counter revCounter;

    /* Adjacent duplicate read dedup flag, with its own adaptive probability. */
    Counter dupCounter;
    /* Previous record's quality string (storage order), kept for dedup. */
    std::vector<uint8_t> prevQual;
    uint32_t prevQualLen = 0;

    Fcv2Cfg cfg;             /* normalized context-parameter tiers */

    /* Whether the fixed alphabet contains the given quality byte; values outside
       it cannot be coded. */
    bool coversByte(uint8_t b) const { return symbolOf[b] >= 0; }

    int alphaSize;
    int symbolOf[256];      /* raw byte -> internal symbol number */
    int byteOf[TREE_CAP];   /* internal symbol number -> raw byte */
    uint16_t quantFreq[TREE_CAP]; /* quantized frequencies, written to the stream so the decoder rebuilds the same tree */
    bool encodeStarted;
    bool flushed;
    bool modelLoaded;

    explicit fcv2_impl(coder_io* ioPtr, const Fcv2Cfg& cfgIn)
        : io(ioPtr), revCounter(COUNTER_INIT),
          dupCounter(COUNTER_INIT), cfg(normalizeCfg(cfgIn)), alphaSize(0),
          encodeStarted(false), flushed(false), modelLoaded(false)
    {
        initTables();
        memset(symbolOf, -1, sizeof(symbolOf));
        memset(byteOf, 0, sizeof(byteOf));
        memset(quantFreq, 0, sizeof(quantFreq));
    }

    /* Combine each model's log-odds prediction with the weights, then convert
       back to a probability. */
    inline int predict(const size_t* base, int node, int wBucket,
                       Counter** slot, int* stretched, int*& wp)
    {
        wp = cm.weightAt(wBucket, node);
        long dot = 0;
        for (int i = 0; i < FCV2_MODEL_COUNT; i++) {
            slot[i] = cm.counterAt(i, base[i], node);
            stretched[i] = g_stretch[counterProb(*slot[i])];
            dot += (long)wp[i] * stretched[i];
        }
        /* squash() already clamps to [1, 4095]. */
        return squash((int)(dot >> 16));
    }

    /* Updates the weights (gradient descent on the coding-cost loss) and each
       model's counter. */
    inline void update(Counter** slot, const int* stretched, int* wp, int p, int isZero)
    {
        int err = (isZero ? 4095 : 0) - p;
        for (int i = 0; i < FCV2_MODEL_COUNT; i++) {
            wp[i] += (err * stretched[i]) >> WEIGHT_LR_SHIFT;
            counterUpdate(*slot[i], isZero);
        }
    }

    /* Converts an in-record position to a sequencing-cycle index (flipped on the
       reverse strand). */
    inline int cycleOf(uint32_t i, uint32_t len, bool rev) const
    {
        uint32_t c = rev ? (len - 1 - i) : i;
        return (c >= (uint32_t)cfg.cycleMax) ? (cfg.cycleMax - 1) : (int)c;
    }

    bool exportModel(std::vector<uint8_t>& out) const
    {
        if (alphaSize <= 0 || alphaSize > TREE_CAP || cm.alphaSize != alphaSize) {
            return false;
        }

        std::vector<uint8_t> snapshot;
        try {
            snapshot.reserve(modelBlobSize(alphaSize, cfg));
            snapshot.insert(snapshot.end(), MODEL_MAGIC, MODEL_MAGIC + sizeof(MODEL_MAGIC));
            appendU16(snapshot, MODEL_FORMAT_VERSION);
            appendU16(snapshot, (uint16_t)cfg.cycleMax);
            appendU16(snapshot, (uint16_t)cfg.cycleBucket);
            appendU16(snapshot, (uint16_t)cfg.prevShift);
            appendU16(snapshot, (uint16_t)(cfg.useDedup ? 1 : 0));
            appendU16(snapshot, (uint16_t)TREE_CAP);
            appendU16(snapshot, (uint16_t)alphaSize);
            appendU32(snapshot, (uint32_t)cm.m0.size());
            appendU32(snapshot, (uint32_t)cm.m1.size());
            appendU32(snapshot, (uint32_t)cm.m2.size());
            appendU32(snapshot, (uint32_t)cm.m3.size());
            appendU32(snapshot, (uint32_t)cm.m4.size());
            appendU32(snapshot, (uint32_t)cm.weight.size());
            for (int s = 0; s < alphaSize; s++) appendU8(snapshot, (uint8_t)byteOf[s]);
            for (int s = 0; s < alphaSize; s++) appendU16(snapshot, quantFreq[s]);
            appendU16(snapshot, revCounter);
            appendU16(snapshot, dupCounter);
            appendCounters(snapshot, cm.m0);
            appendCounters(snapshot, cm.m1);
            appendCounters(snapshot, cm.m2);
            appendCounters(snapshot, cm.m3);
            appendCounters(snapshot, cm.m4);
            for (size_t i = 0; i < cm.weight.size(); i++) appendI32(snapshot, cm.weight[i]);
        } catch (...) {
            return false;
        }
        out.swap(snapshot);
        return true;
    }

    bool loadModel(const std::vector<uint8_t>& blob)
    {
        ModelBlobReader reader(blob);
        for (size_t i = 0; i < sizeof(MODEL_MAGIC); i++) {
            uint8_t byte = 0;
            if (!reader.readU8(byte) || byte != MODEL_MAGIC[i]) return false;
        }

        uint16_t version = 0, cycleMax = 0, cycleBucket = 0;
        uint16_t prevShift = 0, useDedup = 0, treeCap = 0, storedAlpha = 0;
        uint32_t m0Length = 0, m1Length = 0, m2Length = 0, m3Length = 0, m4Length = 0;
        uint32_t weightLength = 0;
        if (!reader.readU16(version) || !reader.readU16(cycleMax) ||
            !reader.readU16(cycleBucket) || !reader.readU16(prevShift) ||
            !reader.readU16(useDedup) || !reader.readU16(treeCap) ||
            !reader.readU16(storedAlpha) ||
            !reader.readU32(m0Length) || !reader.readU32(m1Length) ||
            !reader.readU32(m2Length) || !reader.readU32(m3Length) ||
            !reader.readU32(m4Length) || !reader.readU32(weightLength)) {
            return false;
        }
        if (version != MODEL_FORMAT_VERSION || treeCap != TREE_CAP ||
            storedAlpha == 0 || storedAlpha > TREE_CAP) {
            return false;
        }
        if (useDedup > 1 || prevShift > CFG_PREV_SHIFT_MAX) {
            return false;
        }

        Fcv2Cfg storedCfg = { (int)cycleMax, (int)cycleBucket, (int)prevShift, useDedup == 1 };
        storedCfg = normalizeCfg(storedCfg);
        /* Adopt the prior's own tiers as the effective tiers; a later tier
           mismatch with the stream header is rejected in begin_decode. */
        cfg = storedCfg;

        ContextModel restored;
        try {
            restored.init((int)storedAlpha, storedCfg);
        } catch (...) {
            return false;
        }
        if (m0Length != restored.m0.size() || m1Length != restored.m1.size() ||
            m2Length != restored.m2.size() || m3Length != restored.m3.size() ||
            m4Length != restored.m4.size() || weightLength != restored.weight.size()) {
            return false;
        }

        int restoredByteOf[TREE_CAP];
        int restoredSymbolOf[256];
        uint16_t restoredQuantFreq[TREE_CAP];
        memset(restoredByteOf, 0, sizeof(restoredByteOf));
        memset(restoredSymbolOf, -1, sizeof(restoredSymbolOf));
        memset(restoredQuantFreq, 0, sizeof(restoredQuantFreq));
        for (int s = 0; s < storedAlpha; s++) {
            uint8_t byte = 0;
            if (!reader.readU8(byte) || restoredSymbolOf[byte] >= 0) return false;
            restoredByteOf[s] = byte;
            restoredSymbolOf[byte] = s;
        }
        uint32_t frequency[TREE_CAP];
        for (int s = 0; s < storedAlpha; s++) {
            if (!reader.readU16(restoredQuantFreq[s]) || restoredQuantFreq[s] == 0) return false;
            frequency[s] = restoredQuantFreq[s];
        }

        Counter restoredRevCounter = 0;
        Counter restoredDupCounter = 0;
        if (!reader.readU16(restoredRevCounter) ||
            !reader.readU16(restoredDupCounter) ||
            !readCounters(reader, restored.m0) || !readCounters(reader, restored.m1) ||
            !readCounters(reader, restored.m2) || !readCounters(reader, restored.m3) ||
            !readCounters(reader, restored.m4) ||
            !readWeights(reader, restored.weight) || !reader.atEnd()) {
            return false;
        }

        HuffTree restoredTree;
        restoredTree.build(frequency, (int)storedAlpha);
        cm = std::move(restored);
        cm.rebindTables();   /* see the note on ContextModel::table */
        tree = restoredTree;
        revCounter = restoredRevCounter;
        dupCounter = restoredDupCounter;
        prevQualLen = 0;   /* the prior carries only model state, not the previous quality string */
        prevQual.clear();
        alphaSize = (int)storedAlpha;
        memcpy(byteOf, restoredByteOf, sizeof(byteOf));
        memcpy(symbolOf, restoredSymbolOf, sizeof(symbolOf));
        memcpy(quantFreq, restoredQuantFreq, sizeof(quantFreq));
        modelLoaded = true;
        return true;
    }

private:
    static void appendU8(std::vector<uint8_t>& out, uint8_t value)
    {
        out.push_back(value);
    }

    static void appendCounters(std::vector<uint8_t>& out, const std::vector<Counter>& values)
    {
        for (size_t i = 0; i < values.size(); i++) appendU16(out, values[i]);
    }

    static bool readCounters(ModelBlobReader& reader, std::vector<Counter>& values)
    {
        for (size_t i = 0; i < values.size(); i++) {
            if (!reader.readU16(values[i])) return false;
        }
        return true;
    }

    static bool readWeights(ModelBlobReader& reader, std::vector<int>& values)
    {
        for (size_t i = 0; i < values.size(); i++) {
            uint32_t bits = 0;
            if (!reader.readU32(bits)) return false;
            values[i] = bits <= (uint32_t)std::numeric_limits<int32_t>::max()
                ? (int32_t)bits
                : (int32_t)((int64_t)bits - ((int64_t)1 << 32));
        }
        return true;
    }

    static size_t modelBlobSize(int alpha, const Fcv2Cfg& cfgIn)
    {
        Fcv2Cfg cfg = normalizeCfg(cfgIn);
        size_t prevStates = (size_t)alpha + 1;
        size_t prevQStates = (size_t)(alpha >> cfg.prevShift) + 2;
        size_t counters = (size_t)TREE_CAP +
            (size_t)cfg.cycleMax * TREE_CAP +
            prevStates * cfg.cycleMax * TREE_CAP +
            prevQStates * prevQStates * cfg.cycleBucket * 2 * TREE_CAP +
            (size_t)BASE_STATES * cfg.cycleMax * TREE_CAP;
        /* header: magic + 7 x u16 + 6 x u32 */
        return sizeof(MODEL_MAGIC) + 7 * 2 + 6 * 4 + (size_t)alpha +
            (size_t)alpha * 2 + 2 + 2 + counters * 2 +
            (size_t)cfg.cycleBucket * TREE_CAP * FCV2_MODEL_COUNT * 4;
    }
};

coder_fcv2::coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable)
    : coder_fcv2(io, freqTable, Fcv2Cfg())
{
}

coder_fcv2::coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable, const Fcv2Cfg& cfg)
    : impl(new fcv2_impl(io, cfg))
{
    /* Only quality values that actually occurred are admitted; a smaller
       alphabet gives a shallower tree. */
    int alpha = 0;
    for (size_t b = 0; b < freqTable.size() && b < 256; b++) {
        if (freqTable[b] == 0) {
            continue;
        }
        if (alpha >= TREE_CAP) {
            /* Alphabet exceeds the tree capacity; leave alphaSize at 0 to signal another coder. */
            alpha = 0;
            break;
        }
        impl->symbolOf[b] = alpha;
        impl->byteOf[alpha] = (int)b;
        alpha++;
    }
    impl->alphaSize = alpha;
    if (alpha > 0) {
        /* Quantize frequencies to 1..65535 so encoder and decoder build the
           same tree from the same values. */
        uint32_t maxFreq = 1;
        for (int s = 0; s < alpha; s++) {
            uint32_t f = freqTable[impl->byteOf[s]];
            if (f > maxFreq) maxFreq = f;
        }
        uint32_t ordered[TREE_CAP];
        for (int s = 0; s < alpha; s++) {
            uint64_t scaled = (uint64_t)freqTable[impl->byteOf[s]] * 65535ULL / maxFreq;
            if (scaled < 1) scaled = 1;
            impl->quantFreq[s] = (uint16_t)scaled;
            ordered[s] = impl->quantFreq[s];
        }
        impl->tree.build(ordered, alpha);
        impl->cm.init(alpha, impl->cfg);
    }
}

coder_fcv2::coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable,
                       const std::vector<uint8_t>& modelBlob, bool* modelLoaded)
    : coder_fcv2(io, freqTable, Fcv2Cfg())
{
    bool loaded = impl->loadModel(modelBlob);
    if (modelLoaded != nullptr) {
        *modelLoaded = loaded;
    }
}

coder_fcv2::coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable, const Fcv2Cfg& cfg,
                       const std::vector<uint8_t>& modelBlob, bool* modelLoaded)
    : coder_fcv2(io, freqTable, cfg)
{
    bool loaded = impl->loadModel(modelBlob);
    if (modelLoaded != nullptr) {
        *modelLoaded = loaded;
    }
}

coder_fcv2::~coder_fcv2() = default;

bool coder_fcv2::coversByte(uint8_t b) const
{
    return impl->coversByte(b);
}

bool coder_fcv2::export_model(std::vector<uint8_t>& out) const
{
    return impl->exportModel(out);
}

void coder_fcv2::encode_record(const uint8_t* qual, uint32_t len, bool rev,
                               const uint8_t* seq, uint32_t seqLen)
{
    fcv2_impl* d = impl.get();
    if (d->alphaSize == 0 || qual == nullptr || len == 0) {
        return;
    }
    if (!d->encodeStarted) {
        d->rc.InitEncoder(d->io->data, d->io->data_capacity);
        /* Format version first, then the parameters. */
        d->rc.EncodeByte(FCV2_STREAM_VERSION);
        d->rc.EncodeByte((unsigned)d->alphaSize);
        /* Context parameter tiers, so the decoder can align its context indices. */
        d->rc.EncodeByte((unsigned)d->cfg.cycleMax);
        d->rc.EncodeByte((unsigned)d->cfg.cycleBucket);
        d->rc.EncodeByte((unsigned)d->cfg.prevShift);
        d->rc.EncodeByte((unsigned)(d->cfg.useDedup ? 1 : 0));
        /* Alphabet plus quantized frequencies, so the decoder rebuilds the same tree. */
        for (int s = 0; s < d->alphaSize; s++) {
            d->rc.EncodeByte((unsigned)d->byteOf[s]);
            d->rc.EncodeByte((unsigned)(d->quantFreq[s] >> 8));
            d->rc.EncodeByte((unsigned)(d->quantFreq[s] & 0xFF));
        }
        d->encodeStarted = true;
        d->io->m = coder_io::MENC;
        d->io->appen_magic("coder_fcv2");
    }

    int rv = rev ? 1 : 0;
    {
        int p0 = counterProb(d->revCounter);
        if (p0 < 1)    p0 = 1;
        if (p0 > 4095) p0 = 4095;
        if (rv) {
            d->rc.EncodeBit1<12>(p0);
        } else {
            d->rc.EncodeBit0<12>(p0);
        }
        counterUpdate(d->revCounter, !rv);
    }

    /* Adjacent duplicate read dedup: 1 bit when identical, else the whole string. */
    if (d->cfg.useDedup) {
        bool isDup = (len == d->prevQualLen) &&
                     (d->prevQualLen == 0 || memcmp(qual, d->prevQual.data(), len) == 0);
        int p0 = counterProb(d->dupCounter);
        if (p0 < 1)    p0 = 1;
        if (p0 > 4095) p0 = 4095;
        if (isDup) {
            d->rc.EncodeBit1<12>(p0);
        } else {
            d->rc.EncodeBit0<12>(p0);
        }
        counterUpdate(d->dupCounter, !isDup);
        if (isDup) {
            d->prevQualLen = len;
            if (d->prevQual.size() < (size_t)len) {
                d->prevQual.resize(len);
            }
            memcpy(d->prevQual.data(), qual, len);
            return;
        }
    }

    int prev = d->alphaSize;
    int prev2 = d->alphaSize;

    for (uint32_t i = 0; i < len; i++) {
        int sym = d->symbolOf[qual[i]];
        if (sym < 0) {
            /* A value outside the fixed alphabet cannot be coded; fail the
               stream so the caller falls back to another codec. */
            d->io->set_err(coder_io::IO_UNCODABLE);
            return;
        }
        /* Base of the cycle that produced this quality value; unknown when seq
           is missing or too short. */
        uint32_t mapped = rev ? (len - 1 - i) : i;
        int cyc = d->cycleOf(i, len, rev);
        int baseSym = (seq != nullptr && mapped < seqLen) ? g_baseSym[seq[mapped]] : (BASE_STATES - 1);
        size_t base[FCV2_MODEL_COUNT];
        int wBucket = 0;
        d->cm.slotBases(prev, prev2, cyc, rv, baseSym, base, wBucket);

        int pathLen = d->tree.pathLen[sym];
        for (int step = 0; step < pathLen; step++) {
            int node = d->tree.pathNode[sym][step];
            int bit  = d->tree.pathBit[sym][step];
            Counter* slot[FCV2_MODEL_COUNT];
            int stretched[FCV2_MODEL_COUNT];
            int* wp = nullptr;
            int p0 = d->predict(base, node, wBucket, slot, stretched, wp);
            if (bit) {
                d->rc.EncodeBit1<12>(p0);
            } else {
                d->rc.EncodeBit0<12>(p0);
            }
            d->update(slot, stretched, wp, p0, !bit);
        }
        prev2 = prev;
        prev = sym;
    }

    /* Non-duplicate record: store this record's quality string for the next
       record's dedup comparison. */
    d->prevQualLen = len;
    if (d->prevQual.size() < (size_t)len) {
        d->prevQual.resize(len);
    }
    memcpy(d->prevQual.data(), qual, len);
}

int32_t coder_fcv2::encode_flush()
{
    fcv2_impl* d = impl.get();
    if (!d->encodeStarted || d->flushed) {
        return d->io->data_len;
    }
    d->flushed = true;
    d->io->data_len = d->rc.FinishEncoder();
    return d->io->data_len;
}

/* Reads the stream header, rewinding the range decoder first. Anything that
   does not fit this build's layout is refused as an unsupported version. */
static bool readStreamHeader(FcRangeCoder& rc, const uint8_t* data, int& alpha, Fcv2Cfg& cfg)
{
    rc.InitDecoder(data);
    if ((int)rc.DecodeByte() != (int)FCV2_STREAM_VERSION) {
        return false;
    }
    alpha = (int)rc.DecodeByte();
    if (alpha <= 0 || alpha > TREE_CAP) {
        return false;
    }
    cfg.cycleMax = (int)rc.DecodeByte();
    cfg.cycleBucket = (int)rc.DecodeByte();
    cfg.prevShift = (int)rc.DecodeByte();
    const int useDedup = (int)rc.DecodeByte();
    if (cfg.prevShift > CFG_PREV_SHIFT_MAX || cfg.cycleBucket > CFG_CYCLE_BUCKET_MAX ||
        useDedup > 1) {
        return false;   /* a parameter lies outside its legal range */
    }
    cfg.useDedup = (useDedup != 0);
    return true;
}

int32_t coder_fcv2::begin_decode()
{
    fcv2_impl* d = impl.get();
    /* A stream that does not fit this build's layout is refused as unsupported. */
    int alpha = 0;
    Fcv2Cfg streamCfg;
    if (!readStreamHeader(d->rc, d->io->data, alpha, streamCfg)) {
        return coder_ns::CODER_ERR_UNSUPPORTED_VERSION;
    }
    streamCfg = normalizeCfg(streamCfg);
    /* A prior whose tiers disagree with the header can only be rejected. */
    if (d->modelLoaded && d->cfg != streamCfg) {
        return -1;
    }
    d->cfg = streamCfg;

    int streamByteOf[TREE_CAP];
    uint16_t streamQuantFreq[TREE_CAP];
    memset(streamByteOf, 0, sizeof(streamByteOf));
    memset(streamQuantFreq, 0, sizeof(streamQuantFreq));
    uint32_t freq[TREE_CAP];
    for (int s = 0; s < alpha; s++) {
        int b = (int)d->rc.DecodeByte();
        int hi = (int)d->rc.DecodeByte();
        int lo = (int)d->rc.DecodeByte();
        streamByteOf[s] = b;
        freq[s] = (uint32_t)((hi << 8) | lo);
        streamQuantFreq[s] = (uint16_t)freq[s];
    }

    /* Prior counters are indexed by tree node, so every header item must match
       the loaded prior. */
    if (d->modelLoaded) {
        if (d->alphaSize != alpha) return -1;
        for (int s = 0; s < alpha; s++) {
            if (d->byteOf[s] != streamByteOf[s] || d->quantFreq[s] != streamQuantFreq[s]) {
                return -1;
            }
        }
        return 0;
    }

    memset(d->symbolOf, -1, sizeof(d->symbolOf));
    for (int s = 0; s < alpha; s++) {
        d->byteOf[s] = streamByteOf[s];
        d->symbolOf[streamByteOf[s]] = s;
        d->quantFreq[s] = streamQuantFreq[s];
    }
    d->alphaSize = alpha;
    /* Rebuild from the quantized frequencies carried in the stream, so the
       tree is identical to the encoder's. */
    d->tree.build(freq, alpha);
    d->cm.init(alpha, d->cfg);
    return 0;
}

int32_t coder_fcv2::decode_record(uint8_t* dst, uint32_t len,
                                  const uint8_t* seq, uint32_t seqLen)
{
    fcv2_impl* d = impl.get();
    if (d->alphaSize == 0 || dst == nullptr || len == 0) {
        return 0;
    }

    /* The strand is carried in the stream, see the note on fcv2_impl::revCounter. */
    int p0 = counterProb(d->revCounter);
    if (p0 < 1)    p0 = 1;
    if (p0 > 4095) p0 = 4095;
    int rv = d->rc.DecodeBit<12>(p0);
    counterUpdate(d->revCounter, !rv);
    bool rev = (rv != 0);

    /* Adjacent duplicate read dedup: when dup=1, copy the quality string
       decoded for the previous record directly and skip this record. */
    if (d->cfg.useDedup) {
        int pDup = counterProb(d->dupCounter);
        if (pDup < 1)    pDup = 1;
        if (pDup > 4095) pDup = 4095;
        int isDup = d->rc.DecodeBit<12>(pDup);
        counterUpdate(d->dupCounter, !isDup);
        if (isDup) {
            if (d->prevQualLen < len) {
                /* Corrupt stream: previous string too short; fill with '!' to stay in bounds. */
                memset(dst, '!', len);
            } else {
                memcpy(dst, d->prevQual.data(), len);
            }
            d->prevQualLen = len;
            if (d->prevQual.size() < (size_t)len) {
                d->prevQual.resize(len);
            }
            memcpy(d->prevQual.data(), dst, len);
            return (int32_t)len;
        }
    }

    int prev = d->alphaSize;
    int prev2 = d->alphaSize;

    for (uint32_t i = 0; i < len; i++) {
        uint32_t mapped = rev ? (len - 1 - i) : i;
        int cyc = d->cycleOf(i, len, rev);
        int baseSym = (seq != nullptr && mapped < seqLen) ? g_baseSym[seq[mapped]] : (BASE_STATES - 1);
        size_t base[FCV2_MODEL_COUNT];
        int wBucket = 0;
        d->cm.slotBases(prev, prev2, cyc, rv, baseSym, base, wBucket);

        int cur = d->tree.root;
        while (cur >= d->alphaSize) {
            int node = cur - d->alphaSize;
            Counter* slot[FCV2_MODEL_COUNT];
            int stretched[FCV2_MODEL_COUNT];
            int* wp = nullptr;
            int p0 = d->predict(base, node, wBucket, slot, stretched, wp);
            int bit = d->rc.DecodeBit<12>(p0);
            d->update(slot, stretched, wp, p0, !bit);
            cur = d->tree.child[node][bit];
        }
        dst[i] = (uint8_t)d->byteOf[cur];
        prev2 = prev;
        prev = cur;
    }
    d->prevQualLen = len;
    if (d->prevQual.size() < (size_t)len) {
        d->prevQual.resize(len);
    }
    memcpy(d->prevQual.data(), dst, len);
    return (int32_t)len;
}
