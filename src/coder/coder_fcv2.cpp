/*
 * coder_fcv2.cpp - Implementation of the quality-value context-mixing coder
 *
 * The entire implementation lives in this file; the header only exposes the
 * per-record call interface. See the note in coder_fcv2.h for the reason:
 * this file needs the binary range coder from fc/rangecoder.h, which conflicts
 * with the same-named class in clr.h, and the callers also need clr.h.
 *
 * Algorithm overview:
 *   A quality-value symbol is first split by a Huffman tree into a chain of
 *   binary decisions. At each decision point several probability models of
 *   differing granularity each predict, and the predictions are summed in the
 *   log-odds domain with weights before being handed to the range coder. The
 *   weights are updated by gradient descent; the loss is exactly the number of
 *   bits emitted, so the weights directly optimize the compression ratio.
 */

#include "coder_fcv2.h"

#include <string.h>
#include <math.h>

#include <limits>

#include "coder.h"
#include "fc/rangecoder.h"

namespace {

/* Upper bound on the number of Huffman tree nodes; also limits the alphabet size. */
const int TREE_CAP = 64;

/* Weight-update step size as a shift, equivalent to a learning rate of 1/4096. */
const int WEIGHT_LR_SHIFT = 12;

/* Number of models participating in the mix. m0..m5 are the six existing tiers;
   m6 is the read average-quality tier added by strategy 4. */
const int MODEL_COUNT = 7;

/* Number of read average-quality bins for m6 (strategy 4). */
const int QA_BINS = 4;

/* Legal range of the context parameter tiers. When the decoder reads back the
   stream header it normalizes by the same rule, so that a corrupted stream
   cannot cause an oversized model array to be allocated. Array sizes grow with
   cfg, so an upper bound is required. */
const int CFG_CYCLE_MAX_MIN = 32;
const int CFG_CYCLE_MAX_MAX = 128;
const int CFG_CYCLE_BUCKET_MAX = 32;
const int CFG_DELTA_MAX_MAX = 256;
const int CFG_DELTA_BUCKET_MAX = 16;
const int CFG_PREV_SHIFT_MAX = 3;

/*
 * Normalizes an externally supplied cfg to the tiers actually used by the
 * coder. Both sides (encoder and decoder) must call this; if either side
 * diverges from it, encode/decode will misalign on the context indices.
 */
Fcv2Cfg normalizeCfg(Fcv2Cfg cfg)
{
    if (cfg.cycleMax < CFG_CYCLE_MAX_MIN) cfg.cycleMax = CFG_CYCLE_MAX_MIN;
    if (cfg.cycleMax > CFG_CYCLE_MAX_MAX) cfg.cycleMax = CFG_CYCLE_MAX_MAX;
    if (cfg.cycleBucket < 1) cfg.cycleBucket = 1;
    if (cfg.cycleBucket > CFG_CYCLE_BUCKET_MAX) cfg.cycleBucket = CFG_CYCLE_BUCKET_MAX;
    if (cfg.deltaMax < 8) cfg.deltaMax = 8;
    if (cfg.deltaMax > CFG_DELTA_MAX_MAX) cfg.deltaMax = CFG_DELTA_MAX_MAX;
    if (cfg.deltaBucket < 1) cfg.deltaBucket = 1;
    if (cfg.deltaBucket > CFG_DELTA_BUCKET_MAX) cfg.deltaBucket = CFG_DELTA_BUCKET_MAX;
    if (cfg.prevShift < 0) cfg.prevShift = 0;
    if (cfg.prevShift > CFG_PREV_SHIFT_MAX) cfg.prevShift = CFG_PREV_SHIFT_MAX;
    if (!cfg.useDelta) {
        /* When the delta context is disabled the bin count is normalized to 1,
           the model degenerates to "position + tree node", and the array layout
           is unchanged. */
        cfg.deltaBucket = 1;
    }
    return cfg;
}

/*
 * Computes the average-quality bin of a read (strategy 4, fqzcomp's do_qa).
 * The mean is independent of storage order (the sum runs over the whole quality
 * string); the phred mean is divided into 4 bins. Only the encoder uses it to
 * write the bin into the stream; the decoder reads it back directly, so these
 * thresholds only affect the compression ratio, not both-side consistency.
 */
inline int qaBinOf(const uint8_t* qual, uint32_t len)
{
    uint64_t sum = 0;
    for (uint32_t i = 0; i < len; i++) {
        sum += qual[i];
    }
    int avg = (int)(sum / len);   /* raw byte value, about 33..74 (phred 0..41) */
    if (avg < 63) return 0;       /* phred < 30 */
    if (avg < 68) return 1;       /* 30..34 */
    if (avg < 71) return 2;       /* 35..37 */
    return 3;                     /* >= 38 */
}

/* Number of distinct base-context states: ACGTN + unknown (seq empty or a
   non-base character). */
const int BASE_STATES = 6;

/* Raw byte -> base symbol (ACGTN -> 0..4), everything else maps to 5
   (unknown). Filled in initTables. */
uint8_t g_baseSym[256];

/*
 * The model snapshot is deliberately kept separate from the stream: the stream
 * format must not change because of the prior-model feature. The fixed header
 * below lets the snapshot be validated against the current array layout before
 * it is read; if any size constant changes, the node indices in an old snapshot
 * are no longer well-defined and must be rejected, rather than continuing to
 * encode with probabilities that look usable but are actually wrong.
 */
const uint8_t MODEL_MAGIC[] = { 'f', 'c', 'v', '2', 'p', 'r', 'i', 'o', 'r' };
const uint16_t MODEL_FORMAT_VERSION = 7;  /* v7: read average-quality tier context m6 (useQa + QA_BINS) */

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

/*
 * The byte reader advances its position only when enough input remains, so a
 * truncated snapshot cannot cause an out-of-bounds read, nor can it expose a
 * half-restored model to the caller. All numbers are parsed byte-by-byte in
 * little-endian order, keeping host endianness and struct padding out of the
 * snapshot format.
 */
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

/*
 * A probability and a count are packed into one uint16_t: the low 12 bits are
 * the probability (0..4095) and the high 4 bits are the update count (0..15).
 *
 * The count is capped at 15 rather than something larger because the adaptive
 * step-size table becomes constant once the count reaches 15, so recording
 * anything beyond that is pointless. This conveniently frees 4 bits to share a
 * single 16-bit integer with the probability, halving the model's memory.
 */
typedef uint16_t Counter;

inline int counterProb(Counter c)  { return c & 0x0FFF; }
inline int counterCount(Counter c) { return c >> 12; }

const Counter COUNTER_INIT = 2048;   /* probability 0.5, count 0 */

/* Lookup tables converting between log-odds and probability, 12-bit fixed
   point with a scale factor of 256. */
short  g_stretch[4096];
uint16_t g_squash[4096];

/*
 * Adaptive step-size table.
 *
 * A freshly created counter has few updates and a large step, so a few updates
 * suffice to move it from the initial 0.5 toward the true value; as the update
 * count grows the step shrinks and the estimate gradually stabilizes. This
 * targets sparse contexts directly - the finer the context, the fewer times
 * each counter is accessed, and it sits exactly in the phase of fast convergence
 * with a large step.
 */
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

/*
 * Huffman tree.
 *
 * The leaves are quality-value symbols; the path from root to leaf is a chain
 * of binary decisions. Building it by frequency gives high-frequency symbols
 * shorter paths, reducing the number of binary codings per symbol. This only
 * affects speed, not the compression ratio: the cost of arithmetic coding is
 * the sum of -log2(p) over the steps, and splitting a symbol into several
 * conditional decisions leaves the total cost unchanged.
 */
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
 * Seven context models of increasing granularity.
 *
 *   m0  only the tree node; the coarsest, always has enough samples, acts as a
 *       floor/baseline
 *   m1  position + tree node
 *   m2  previous symbol + position + tree node
 *   m3  previous two symbols + position bin + strand + tree node
 *   m4  current base + position + tree node
 *   m5  in-read quality-transition count bin + position + tree node
 *       (strategy 1)
 *   m6  read average-quality bin + position + tree node (strategy 4)
 *
 * The base used by m4 is the base of "the sequencing cycle that produced this
 * quality value": position i on the forward strand, len-1-i in storage order on
 * the reverse strand. Quality values correlate strongly with the base
 * (mismatch positions and low-complexity regions are systematically lower in
 * quality), and this context dimension is largely orthogonal to the others
 * (previous symbols, strand), so it does not overlap m3.
 *
 * m5's transition count is another description of quality monotonically
 * changing along a read (later sequencing cycles tend to have lower quality):
 * "how many times the quality has changed up to the current position" reflects
 * whether this read's quality is clean, does not overlap the local previous
 * symbols (m2/m3) or the position (m1), and accumulates only within a read,
 * resetting per record. fqzcomp uses the same dimension (state->delta) in its
 * context; here it is made an independent model handed to the mixer, and the
 * weights decide how much it is worth.
 *
 * m6 is the whole-read average-quality bin (strategy 4, fqzcomp's do_qa):
 * average quality varies widely across reads (empirically a long tail from
 * Q10 to Q40), so whether an entire read sits in the "high-quality bin" or the
 * "low-quality bin" is a strong prior on the symbol distribution. The bin is
 * computed by the encoder and written into the stream (2 bits per record); the
 * decoder reads it back directly, keeping both sides consistent.
 *
 * The sizes of the arrays are determined by Fcv2Cfg, see init().
 */
struct ContextModel {
    int alphaSize;
    int prevStates;          /* number of symbol values plus one; the extra slot represents "no predecessor" */
    int prevQStates;         /* number of quantized predecessor values for m3 (>> prevShift) */
    int cycleMax;
    int cycleBucket;
    int deltaMax;
    int deltaBucket;
    int prevShift;
    std::vector<Counter> m0, m1, m2, m3, m4, m5, m6;
    std::vector<int> weight; /* mixing weights, organized by position bin x tree node x model */

    void init(int alpha, const Fcv2Cfg& cfg)
    {
        alphaSize = alpha;
        prevStates = alpha + 1;
        cycleMax = cfg.cycleMax;
        cycleBucket = cfg.cycleBucket;
        deltaMax = cfg.deltaMax;
        deltaBucket = cfg.deltaBucket;
        prevShift = cfg.prevShift;
        /*
         * Using all symbol values for m3 is too sparse (39x39x16x2x64 = about
         * 3.1M slots on this data, only ~9 accesses per block), so the order-2
         * structure cannot be learned. Empirically, quantizing the previous
         * quality value by >>1 to 22 tiers (quality value ~ symbol id) retains
         * 55% of the order-2 conditional-entropy gain (3.66 vs 3.53 bits) while
         * cutting the slot count to 1/3. The tier count is adjustable via
         * prevShift: shift further right for coarser tiers with large
         * alphabets, or skip quantization for small alphabets.
         */
        prevQStates = (alpha >> prevShift) + 2;
        m0.assign((size_t)TREE_CAP, COUNTER_INIT);
        m1.assign((size_t)cycleMax * TREE_CAP, COUNTER_INIT);
        m2.assign((size_t)prevStates * cycleMax * TREE_CAP, COUNTER_INIT);
        m3.assign((size_t)prevQStates * prevQStates * cycleBucket * 2 * TREE_CAP, COUNTER_INIT);
        m4.assign((size_t)BASE_STATES * cycleMax * TREE_CAP, COUNTER_INIT);
        m5.assign((size_t)deltaBucket * cycleMax * TREE_CAP, COUNTER_INIT);
        m6.assign((size_t)QA_BINS * cycleMax * TREE_CAP, COUNTER_INIT);
        /* Initial weights are 1<<14, i.e. fixed-point 0.25, so the seven models
           together start close to a simple average. */
        weight.assign((size_t)cycleBucket * TREE_CAP * MODEL_COUNT, 1 << 14);
    }

    inline void slotBases(int prev, int prev2, int cyc, int rev, int baseSym,
                          int delta, int qa, size_t* base) const
    {
        int bucket = cyc * cycleBucket / cycleMax;
        if (bucket >= cycleBucket) bucket = cycleBucket - 1;
        int qp = (prev >> prevShift);      /* quantized previous quality value */
        int qp2 = (prev2 >> prevShift);
        if (qp >= prevQStates) qp = prevQStates - 1;
        if (qp2 >= prevQStates) qp2 = prevQStates - 1;
        int dlt = delta >= deltaMax ? deltaBucket - 1 : delta * deltaBucket / deltaMax;
        base[0] = 0;
        base[1] = (size_t)cyc * TREE_CAP;
        base[2] = ((size_t)prev * cycleMax + cyc) * TREE_CAP;
        base[3] = ((((size_t)qp2 * prevQStates + qp) * cycleBucket + bucket) * 2 + rev) * TREE_CAP;
        base[4] = ((size_t)baseSym * cycleMax + cyc) * TREE_CAP;
        base[5] = ((size_t)dlt * cycleMax + cyc) * TREE_CAP;
        base[6] = ((size_t)qa * cycleMax + cyc) * TREE_CAP;
    }

    inline int* weightPtr(int cyc, int node)
    {
        int bucket = cyc * cycleBucket / cycleMax;
        if (bucket >= cycleBucket) bucket = cycleBucket - 1;
        return &weight[((size_t)bucket * TREE_CAP + node) * MODEL_COUNT];
    }

    inline Counter* counterAt(int model, size_t base, int node)
    {
        switch (model) {
        case 0:  return &m0[base + node];
        case 1:  return &m1[base + node];
        case 2:  return &m2[base + node];
        case 3:  return &m3[base + node];
        case 4:  return &m4[base + node];
        case 5:  return &m5[base + node];
        default: return &m6[base + node];
        }
    }
};

} /* namespace */

/*
 * The implementation body. It lives outside the anonymous namespace because
 * the header forward-declares it.
 */
class fcv2_impl {
public:
    coder_io* io;
    ContextModel cm;
    HuffTree tree;
    RangeCoder rc;

    /*
     * Probability model for the strand itself.
     *
     * The strand is written into the stream with each record, rather than being
     * supplied again by the caller at decompression. This makes the coder fully
     * self-contained: the decompression side need not care how far the FLAG
     * field has been decoded in order to decode QUAL, nor maintain a
     * per-record array of strands.
     *
     * The cost is about one bit per record. Empirically the reverse strand is
     * 49.9%, close to maximum entropy, so it compresses poorly: a million
     * records cost about 125 KB, which is 0.139% of 90 MB of quality values.
     * Meanwhile the strand brings a gain of 0.372 percentage points, so after
     * subtracting this cost the net gain is still 0.233 percentage points;
     * worth it.
     *
     * An adaptive counter is still used instead of a fixed 0.5 because forward
     * and reverse strands are not necessarily split evenly in a coordinate-
     * sorted SAM; the model adapts to the actual ratio.
     */
    Counter revCounter;

    /* Adjacent duplicate read dedup (strategy 3): the dup flag has its own
       adaptive probability, maintained record by record. */
    Counter dupCounter;
    /* The previous record's quality string (storage order), used for
       comparison and for the decoder to restore it; kept across records. */
    std::vector<uint8_t> prevQual;
    uint32_t prevQualLen = 0;

    Fcv2Cfg cfg;             /* normalized context-parameter tiers */

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

    /*
     * A single binary prediction: take each model's probability to log-odds,
     * sum them weighted by the weights, then convert back to a probability.
     *
     * The linear combination is done in the log-odds domain rather than as a
     * direct weighted average of probabilities because probability scales are
     * not uniform: the difference between 0.50 and 0.51 and the difference
     * between 0.98 and 0.99 are entirely different in coding cost. Log-odds is
     * an additive evidence strength, so the combination semantics are correct.
     */
    inline int predict(const size_t* base, int node, int cyc,
                       Counter** slot, int* stretched, int*& wp)
    {
        wp = cm.weightPtr(cyc, node);
        long dot = 0;
        for (int i = 0; i < MODEL_COUNT; i++) {
            slot[i] = cm.counterAt(i, base[i], node);
            stretched[i] = g_stretch[counterProb(*slot[i])];
            dot += (long)wp[i] * stretched[i];
        }
        int p = squash((int)(dot >> 16));
        if (p < 1)    p = 1;
        if (p > 4095) p = 4095;
        return p;
    }

    /*
     * Updates the weights and each model's counter.
     *
     * The weight update is gradient descent on the cross-entropy loss. Viewing
     * the mixer as a logistic regression whose inputs are the log-odds given by
     * each model and whose output is the final probability, the gradient of the
     * loss w.r.t. a weight is exactly -(target - prediction) x input. Here the
     * loss is the number of bits actually spent coding this bit, so the weights
     * directly optimize the compression ratio.
     */
    inline void update(Counter** slot, const int* stretched, int* wp, int p, int isZero)
    {
        int err = (isZero ? 4095 : 0) - p;
        for (int i = 0; i < MODEL_COUNT; i++) {
            wp[i] += (err * stretched[i]) >> WEIGHT_LR_SHIFT;
            counterUpdate(*slot[i], isZero);
        }
    }

    /* Converts an in-record position to a sequencing-cycle index; the reverse
       strand must flip it by the record length, see the header note. */
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
            appendU16(snapshot, (uint16_t)cfg.deltaMax);
            appendU16(snapshot, (uint16_t)cfg.deltaBucket);
            appendU16(snapshot, (uint16_t)cfg.prevShift);
            appendU16(snapshot, (uint16_t)(cfg.useDelta ? 1 : 0));
            appendU16(snapshot, (uint16_t)(cfg.useDedup ? 1 : 0));
            appendU16(snapshot, (uint16_t)(cfg.useQa ? 1 : 0));
            appendU16(snapshot, (uint16_t)TREE_CAP);
            appendU16(snapshot, (uint16_t)MODEL_COUNT);
            appendU16(snapshot, (uint16_t)alphaSize);
            appendU32(snapshot, (uint32_t)cm.m0.size());
            appendU32(snapshot, (uint32_t)cm.m1.size());
            appendU32(snapshot, (uint32_t)cm.m2.size());
            appendU32(snapshot, (uint32_t)cm.m3.size());
            appendU32(snapshot, (uint32_t)cm.m4.size());
            appendU32(snapshot, (uint32_t)cm.m5.size());
            appendU32(snapshot, (uint32_t)cm.m6.size());
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
            appendCounters(snapshot, cm.m5);
            appendCounters(snapshot, cm.m6);
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

        uint16_t version = 0, cycleMax = 0, cycleBucket = 0, deltaMax = 0, deltaBucket = 0;
        uint16_t prevShift = 0, useDelta = 0, useDedup = 0, useQa = 0, treeCap = 0, modelCount = 0, storedAlpha = 0;
        uint32_t m0Length = 0, m1Length = 0, m2Length = 0, m3Length = 0, m4Length = 0;
        uint32_t m5Length = 0, m6Length = 0, weightLength = 0;
        if (!reader.readU16(version) || !reader.readU16(cycleMax) ||
            !reader.readU16(cycleBucket) || !reader.readU16(deltaMax) ||
            !reader.readU16(deltaBucket) || !reader.readU16(prevShift) ||
            !reader.readU16(useDelta) || !reader.readU16(useDedup) ||
            !reader.readU16(useQa) ||
            !reader.readU16(treeCap) || !reader.readU16(modelCount) ||
            !reader.readU16(storedAlpha) ||
            !reader.readU32(m0Length) || !reader.readU32(m1Length) ||
            !reader.readU32(m2Length) || !reader.readU32(m3Length) ||
            !reader.readU32(m4Length) || !reader.readU32(m5Length) ||
            !reader.readU32(m6Length) || !reader.readU32(weightLength)) {
            return false;
        }
        if (version != MODEL_FORMAT_VERSION || treeCap != TREE_CAP ||
            modelCount != MODEL_COUNT || storedAlpha == 0 || storedAlpha > TREE_CAP) {
            return false;
        }
        if (useDelta > 1 || useDedup > 1 || useQa > 1 || prevShift > CFG_PREV_SHIFT_MAX) {
            return false;
        }

        Fcv2Cfg storedCfg = { (int)cycleMax, (int)cycleBucket, (int)deltaMax,
                              (int)deltaBucket, (int)prevShift, useDelta == 1,
                              useDedup == 1, useQa == 1 };
        storedCfg = normalizeCfg(storedCfg);
        /*
         * The prior is a cross-block training product, so its counter arrays
         * follow the tier layout used during training. The loader (whether the
         * compression or the decompression side) adopts the tiers carried by
         * the prior itself as the effective tiers, rather than requiring the
         * prior to match the tiers of construction: the prior always shares the
         * tiers of the stream it is used with (training and compression use the
         * same parameter set), so either side can safely treat "the prior
         * decides the tiers". If they disagree with the stream header,
         * begin_decode rejects on the tier mismatch, see its modelLoaded
         * branch.
         */
        cfg = storedCfg;

        ContextModel restored;
        try {
            restored.init((int)storedAlpha, storedCfg);
        } catch (...) {
            return false;
        }
        if (m0Length != restored.m0.size() || m1Length != restored.m1.size() ||
            m2Length != restored.m2.size() || m3Length != restored.m3.size() ||
            m4Length != restored.m4.size() || m5Length != restored.m5.size() ||
            m6Length != restored.m6.size() || weightLength != restored.weight.size()) {
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
            !readCounters(reader, restored.m4) || !readCounters(reader, restored.m5) ||
            !readCounters(reader, restored.m6) ||
            !readWeights(reader, restored.weight) || !reader.atEnd()) {
            return false;
        }

        HuffTree restoredTree;
        restoredTree.build(frequency, (int)storedAlpha);
        cm = std::move(restored);
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
            (size_t)BASE_STATES * cfg.cycleMax * TREE_CAP +
            (size_t)cfg.deltaBucket * cfg.cycleMax * TREE_CAP +
            (size_t)QA_BINS * cfg.cycleMax * TREE_CAP;
        return sizeof(MODEL_MAGIC) + 12 * 2 + 8 * 4 + (size_t)alpha +
            (size_t)alpha * 2 + 2 + 2 + counters * 2 +
            (size_t)cfg.cycleBucket * TREE_CAP * MODEL_COUNT * 4;
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
            /* The alphabet exceeds the tree capacity; the caller should switch
               to another coder. alphaSize is left at 0 so that calls beyond
               supports() are noticed by the upper layer. */
            alpha = 0;
            break;
        }
        impl->symbolOf[b] = alpha;
        impl->byteOf[alpha] = (int)b;
        alpha++;
    }
    impl->alphaSize = alpha;
    if (alpha > 0) {
        /*
         * Frequencies are quantized to 1..65535 before building the tree,
         * rather than using the raw counts directly.
         *
         * This is so that the encoder and decoder see exactly the same input
         * when building the tree: the decoder only gets the quantized values
         * written into the stream, and if the encoder built the tree from raw
         * counts the two trees could diverge at some merge order. Using the
         * quantized values uniformly avoids this problem.
         */
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
        d->rc.EncodeByte((unsigned)d->alphaSize);
        /*
         * The context parameter tiers are written into the stream header. The
         * model array layout is determined by them, and the decoder must read
         * back the same set to align the context indices; this also lets each
         * data block carry its own tier choice. See normalizeCfg(): both sides
         * normalize by the same rule, so any legal value agrees on both sides.
         */
        d->rc.EncodeByte((unsigned)d->cfg.cycleMax);
        d->rc.EncodeByte((unsigned)d->cfg.cycleBucket);
        d->rc.EncodeByte((unsigned)d->cfg.deltaMax);
        d->rc.EncodeByte((unsigned)d->cfg.deltaBucket);
        d->rc.EncodeByte((unsigned)d->cfg.prevShift);
        d->rc.EncodeByte((unsigned)(d->cfg.useDelta ? 1 : 0));
        d->rc.EncodeByte((unsigned)(d->cfg.useDedup ? 1 : 0));
        d->rc.EncodeByte((unsigned)(d->cfg.useQa ? 1 : 0));
        /*
         * The alphabet is written into the stream together with each symbol's
         * quantized frequency.
         *
         * The frequencies must be written out; the decoder must not be left to
         * guess them. The shape of the Huffman tree is entirely determined by
         * the frequencies, and if the trees built on the two sides differ in
         * the slightest way the paths no longer match and the decoded symbol is
         * a different one. After quantization to 16 bits each symbol takes two
         * bytes; a few dozen symbols total under a hundred bytes, negligible
         * relative to a whole block of data.
         */
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

    /*
     * Adjacent duplicate read dedup (strategy 3, fqzcomp's do_dedup): compare
     * byte-by-byte with the previous record; if identical, write 1 bit and skip
     * the whole quality string. When a duplicate is hit, this record's context
     * state is not updated (prev/prev2/delta are all reset at the start of each
     * record), so the next record starts from zero as usual on both sides.
     */
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

    /*
     * read average-quality bin (strategy 4): the bin is written into the stream
     * as the m6 context. Fixed 2 bits/record; the decoder reads back the same
     * bin.
     */
    int qa = 0;
    if (d->cfg.useQa) {
        qa = qaBinOf(qual, len);
        if (qa & 2) {
            d->rc.EncodeBit1<12>(2048);
        } else {
            d->rc.EncodeBit0<12>(2048);
        }
        if (qa & 1) {
            d->rc.EncodeBit1<12>(2048);
        } else {
            d->rc.EncodeBit0<12>(2048);
        }
    }

    int prev = d->alphaSize;
    int prev2 = d->alphaSize;
    int delta = 0;   /* in-read quality-transition count, reset per record */

    for (uint32_t i = 0; i < len; i++) {
        int sym = d->symbolOf[qual[i]];
        if (sym < 0) {
            /* A quality value not seen during the statistics phase cannot be
               coded. The upper layer guarantees this never happens; if it does,
               skipping the byte would break the round trip, and substituting
               the last symbol is likewise unacceptable, so here we keep the
               predecessors unchanged and skip the byte, letting the upper
               layer's round-trip verification expose the problem. */
            continue;
        }
        /*
         * The base context takes the base of "the cycle that produced this
         * quality value": i on the forward strand, len-1-i on the reverse
         * strand, using the same mapping as cycleOf so both sides stay
         * consistent. If the seq length is insufficient (abnormal data) the
         * position falls into the "unknown" bin, and seq is never read
         * out-of-bounds.
         */
        uint32_t mapped = rev ? (len - 1 - i) : i;
        int cyc = d->cycleOf(i, len, rev);
        int baseSym = (seq != nullptr && mapped < seqLen) ? g_baseSym[seq[mapped]] : (BASE_STATES - 1);
        size_t base[MODEL_COUNT];
        d->cm.slotBases(prev, prev2, cyc, rv, baseSym, delta, qa, base);

        int pathLen = d->tree.pathLen[sym];
        for (int step = 0; step < pathLen; step++) {
            int node = d->tree.pathNode[sym][step];
            int bit  = d->tree.pathBit[sym][step];
            Counter* slot[MODEL_COUNT];
            int stretched[MODEL_COUNT];
            int* wp = nullptr;
            int p0 = d->predict(base, node, cyc, slot, stretched, wp);
            if (bit) {
                d->rc.EncodeBit1<12>(p0);
            } else {
                d->rc.EncodeBit0<12>(p0);
            }
            d->update(slot, stretched, wp, p0, !bit);
        }
        /*
         * The transition count is updated only after the symbol is coded, so
         * the m5 context counts the changes "before the current symbol",
         * matching fqzcomp's use of delta. The first symbol has no predecessor
         * and does not count as a transition.
         */
        if (prev != d->alphaSize && sym != prev) {
            delta++;
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

int32_t coder_fcv2::begin_decode()
{
    fcv2_impl* d = impl.get();
    d->rc.InitDecoder(d->io->data);
    int alpha = (int)d->rc.DecodeByte();
    if (alpha <= 0 || alpha > TREE_CAP) {
        return -1;
    }
    /* Read back the context parameter tiers written by the encoder, in the
       exact order encode_record writes them. */
    Fcv2Cfg streamCfg;
    streamCfg.cycleMax = (int)d->rc.DecodeByte();
    streamCfg.cycleBucket = (int)d->rc.DecodeByte();
    streamCfg.deltaMax = (int)d->rc.DecodeByte();
    streamCfg.deltaBucket = (int)d->rc.DecodeByte();
    streamCfg.prevShift = (int)d->rc.DecodeByte();
    streamCfg.useDelta = (d->rc.DecodeByte() != 0);
    streamCfg.useDedup = (d->rc.DecodeByte() != 0);
    streamCfg.useQa = (d->rc.DecodeByte() != 0);
    if (streamCfg.prevShift > CFG_PREV_SHIFT_MAX ||
        streamCfg.cycleBucket > CFG_CYCLE_BUCKET_MAX ||
        streamCfg.deltaBucket > CFG_DELTA_BUCKET_MAX) {
        return -1;   /* corrupted stream: tier exceeds the legal range; refuse to allocate a model from it */
    }
    streamCfg = normalizeCfg(streamCfg);
    /*
     * The prior comes from cross-block training, so its model layout is
     * determined by the tiers used at training time; if it disagrees with the
     * stream header, the counter arrays no longer line up and it can only be
     * rejected (the same mismatch on the compression side makes loadModel
     * return false).
     */
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

    /*
     * The counters in the prior are indexed by tree-node number, so they cannot
     * be used to decode a stream built from a different tree. Even if two
     * alphabets contain the same set of bytes, if the frequencies change and
     * the merge order changes, a node number may now denote a different
     * decision; therefore, after loading a prior, every item in the stream
     * header must be verified, and the loaded learning state must be kept
     * rather than overwritten by calling init again.
     */
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
                /* The previous quality string is missing or too short, a
                   corrupted stream; fall back to '!' to avoid an out-of-bounds
                   read. */
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
    int delta = 0;   /* in-read quality-transition count, maintained in the same order as the encoder */

    /* read average-quality bin: read back the 2 bits written by the encoder,
       as the m6 context. */
    int qa = 0;
    if (d->cfg.useQa) {
        qa = (d->rc.DecodeBit<12>(2048) << 1) | d->rc.DecodeBit<12>(2048);
    }

    for (uint32_t i = 0; i < len; i++) {
        uint32_t mapped = rev ? (len - 1 - i) : i;
        int cyc = d->cycleOf(i, len, rev);
        int baseSym = (seq != nullptr && mapped < seqLen) ? g_baseSym[seq[mapped]] : (BASE_STATES - 1);
        size_t base[MODEL_COUNT];
        d->cm.slotBases(prev, prev2, cyc, rv, baseSym, delta, qa, base);

        int cur = d->tree.root;
        while (cur >= d->alphaSize) {
            int node = cur - d->alphaSize;
            Counter* slot[MODEL_COUNT];
            int stretched[MODEL_COUNT];
            int* wp = nullptr;
            int p0 = d->predict(base, node, cyc, slot, stretched, wp);
            int bit = d->rc.DecodeBit<12>(p0);
            d->update(slot, stretched, wp, p0, !bit);
            cur = d->tree.child[node][bit];
        }
        dst[i] = (uint8_t)d->byteOf[cur];
        if (prev != d->alphaSize && cur != prev) {
            delta++;
        }
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
