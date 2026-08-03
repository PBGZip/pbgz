/*
 * coder_fcv2.cpp - 质量值上下文混合编码器的实现
 *
 * 实现全部集中在本文件，头文件只暴露按记录调用的接口。原因见 coder_fcv2.h 的说明：
 * 本文件要用 fc/rangecoder.h 的二值区间编码器，它与 clr.h 的同名类冲突，而使用方
 * 同时需要 clr.h。
 *
 * 算法梗概：
 *   质量值符号先经哈夫曼树拆成一串二值判断，每个判断点由若干个不同粒度的概率模型
 *   分别预测，预测结果在对数几率域加权求和后交给区间编码器。权重按梯度下降更新，
 *   损失函数就是压缩后的比特数本身，所以权重是在直接优化压缩率。
 */

#include "coder_fcv2.h"

#include <string.h>
#include <math.h>

#include "coder.h"
#include "fc/rangecoder.h"

namespace {

/* 每条记录内的位置（测序循环序号）上限，超出的位置一律归入最后一档。 */
const int CYCLE_MAX = 96;

/* 高阶模型里位置被压缩成多少档。位置本身取值可达 96，直接用会让上下文数量过大。 */
const int CYCLE_BUCKET = 16;

/* 哈夫曼树的节点数上限，同时也限制了字母表大小。 */
const int TREE_CAP = 64;

/* 权重更新的步长移位量，等效学习率 1/4096。 */
const int WEIGHT_LR_SHIFT = 12;

/* 参与混合的模型个数。 */
const int MODEL_COUNT = 4;

/*
 * 概率与计数打包进一个 uint16_t：低 12 位是概率（0..4095），高 4 位是更新次数（0..15）。
 *
 * 计数上限取 15 而不是更大，是因为自适应步长表在计数达到 15 之后就恒定不变了，
 * 再往上记录没有意义。这样正好省出 4 位与概率共用一个 16 位整数，模型内存减半。
 */
typedef uint16_t Counter;

inline int counterProb(Counter c)  { return c & 0x0FFF; }
inline int counterCount(Counter c) { return c >> 12; }

const Counter COUNTER_INIT = 2048;   /* 概率 0.5，计数 0 */

/* 对数几率与概率的互转表，12 位定点，缩放因子 256。 */
short  g_stretch[4096];
uint16_t g_squash[4096];

/*
 * 自适应步长表。
 *
 * 新建的计数器更新次数少，步长大，几次更新就能从初值 0.5 靠近真实值；随着更新次数
 * 增加步长变小，估计逐渐稳定。这直接针对稀疏上下文——上下文越细，每个计数器被访问
 * 的次数越少，正好处在大步长快速收敛的阶段。
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
 * 哈夫曼树。
 *
 * 叶子是质量值符号，从根到叶的路径就是一串二值判断。按频率构造使高频符号路径更短，
 * 减少每符号的二值编码次数。这只影响速度，不影响压缩率：算术编码的代价是各步
 * -log2(p) 之和，把一个符号拆成几步条件判断，总代价不变。
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
            /* 频率加一，保证没出现过的符号也有一条路径，避免编码时越界。 */
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
 * 四个粒度递增的上下文模型。
 *
 *   m0  只看树节点，最粗，任何时候样本都充足，起保底作用
 *   m1  位置 + 树节点
 *   m2  前一个符号 + 位置 + 树节点
 *   m3  前两个符号 + 位置档 + 链方向 + 树节点，最细也最稀疏
 *
 * 早期版本在 m2 和 m3 之间还有一个"前两个符号 + 位置档"的模型，实测它与 m3 的上下文
 * 高度重叠（m3 只是多了链方向一维），信息冗余，去掉后压缩率几乎不变而速度和内存都
 * 改善，因此不再保留。
 */
struct ContextModel {
    int alphaSize;
    int prevStates;          /* 符号取值数加一，多出来的那个表示"没有前驱" */
    std::vector<Counter> m0, m1, m2, m3;
    std::vector<int> weight; /* 混合权重，按 位置档 × 树节点 × 模型 组织 */

    void init(int alpha)
    {
        alphaSize = alpha;
        prevStates = alpha + 1;
        m0.assign((size_t)TREE_CAP, COUNTER_INIT);
        m1.assign((size_t)CYCLE_MAX * TREE_CAP, COUNTER_INIT);
        m2.assign((size_t)prevStates * CYCLE_MAX * TREE_CAP, COUNTER_INIT);
        m3.assign((size_t)prevStates * prevStates * CYCLE_BUCKET * 2 * TREE_CAP, COUNTER_INIT);
        /* 权重初值取 1<<14，即定点 0.25，四个模型合起来接近简单平均。 */
        weight.assign((size_t)CYCLE_BUCKET * TREE_CAP * MODEL_COUNT, 1 << 14);
    }

    inline void slotBases(int prev, int prev2, int cyc, int rev, size_t* base) const
    {
        int bucket = cyc * CYCLE_BUCKET / CYCLE_MAX;
        if (bucket >= CYCLE_BUCKET) bucket = CYCLE_BUCKET - 1;
        base[0] = 0;
        base[1] = (size_t)cyc * TREE_CAP;
        base[2] = ((size_t)prev * CYCLE_MAX + cyc) * TREE_CAP;
        base[3] = ((((size_t)prev2 * prevStates + prev) * CYCLE_BUCKET + bucket) * 2 + rev) * TREE_CAP;
    }

    inline int* weightPtr(int cyc, int node)
    {
        int bucket = cyc * CYCLE_BUCKET / CYCLE_MAX;
        if (bucket >= CYCLE_BUCKET) bucket = CYCLE_BUCKET - 1;
        return &weight[((size_t)bucket * TREE_CAP + node) * MODEL_COUNT];
    }

    inline Counter* counterAt(int model, size_t base, int node)
    {
        switch (model) {
        case 0:  return &m0[base + node];
        case 1:  return &m1[base + node];
        case 2:  return &m2[base + node];
        default: return &m3[base + node];
        }
    }
};

} /* namespace */

/*
 * 实现体。放在匿名命名空间之外是因为头文件里前置声明了它。
 */
class fcv2_impl {
public:
    coder_io* io;
    ContextModel cm;
    HuffTree tree;
    RangeCoder rc;

    int alphaSize;
    int symbolOf[256];      /* 原始字节 -> 内部符号编号 */
    int byteOf[TREE_CAP];   /* 内部符号编号 -> 原始字节 */
    uint16_t quantFreq[TREE_CAP]; /* 量化后的频率，随码流写出供解码端重建同一棵树 */
    bool encodeStarted;
    bool flushed;

    explicit fcv2_impl(coder_io* ioPtr)
        : io(ioPtr), alphaSize(0), encodeStarted(false), flushed(false)
    {
        initTables();
        memset(symbolOf, -1, sizeof(symbolOf));
        memset(byteOf, 0, sizeof(byteOf));
        memset(quantFreq, 0, sizeof(quantFreq));
    }

    /*
     * 一次二值预测：把各模型的概率取对数几率，按权重加权求和，再转回概率。
     *
     * 之所以在对数几率域做线性组合而不是直接对概率加权平均，是因为概率的尺度不均匀：
     * 0.50 和 0.51 之间的差别，与 0.98 和 0.99 之间的差别，在编码代价上完全不是一
     * 回事。对数几率是可加的证据强度，合并语义正确。
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
     * 更新权重和各模型的计数器。
     *
     * 权重更新是对交叉熵损失的梯度下降。把混合器看作一个逻辑回归，输入是各模型给出的
     * 对数几率，输出是最终概率，那么损失对权重的梯度恰好是 -(目标 - 预测) × 输入。
     * 这里的损失函数就是编码这一位实际花掉的比特数，所以权重直接在优化压缩率。
     */
    inline void update(Counter** slot, const int* stretched, int* wp, int p, int isZero)
    {
        int err = (isZero ? 4095 : 0) - p;
        for (int i = 0; i < MODEL_COUNT; i++) {
            wp[i] += (err * stretched[i]) >> WEIGHT_LR_SHIFT;
            counterUpdate(*slot[i], isZero);
        }
    }

    /* 记录内位置转成测序循环序号，反向链需要按记录长度翻转，见头文件说明。 */
    static inline int cycleOf(uint32_t i, uint32_t len, bool rev)
    {
        uint32_t c = rev ? (len - 1 - i) : i;
        return (c >= (uint32_t)CYCLE_MAX) ? (CYCLE_MAX - 1) : (int)c;
    }
};

coder_fcv2::coder_fcv2(coder_io* io, const std::vector<uint32_t>& freqTable)
    : impl(new fcv2_impl(io))
{
    /* 只收录实际出现过的质量值，字母表越小树越浅。 */
    int alpha = 0;
    for (size_t b = 0; b < freqTable.size() && b < 256; b++) {
        if (freqTable[b] == 0) {
            continue;
        }
        if (alpha >= TREE_CAP) {
            /* 字母表超出树容量，调用方应改用其他编码器，这里保持 alphaSize 为 0
               让 supports 之外的调用能被上层察觉。 */
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
         * 频率量化到 1..65535 再建树，而不是直接用原始计数。
         *
         * 这样做是为了让编码端和解码端建树时看到的输入完全一致：解码端只能拿到写进
         * 码流的量化值，如果编码端用原始计数建树，两棵树可能在某个合并顺序上分岔。
         * 统一用量化值就不存在这个问题。
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
        impl->cm.init(alpha);
    }
}

coder_fcv2::~coder_fcv2() = default;

void coder_fcv2::encode_record(const uint8_t* qual, uint32_t len, bool rev)
{
    fcv2_impl* d = impl.get();
    if (d->alphaSize == 0 || qual == nullptr || len == 0) {
        return;
    }
    if (!d->encodeStarted) {
        d->rc.InitEncoder(d->io->data, d->io->data_capacity);
        d->rc.EncodeByte((unsigned)d->alphaSize);
        /*
         * 字母表连同各符号的量化频率一起写进码流。
         *
         * 频率必须写出去，不能让解码端猜。哈夫曼树的形状完全由频率决定，两端建出的
         * 树只要有一点不同，路径就对不上，解出来的就是另一个符号。量化到 16 位后
         * 每个符号两字节，几十个符号总共不到一百字节，相对整块数据可以忽略。
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

    int prev = d->alphaSize;
    int prev2 = d->alphaSize;
    int rv = rev ? 1 : 0;

    for (uint32_t i = 0; i < len; i++) {
        int sym = d->symbolOf[qual[i]];
        if (sym < 0) {
            /* 出现了统计阶段没见过的质量值，无法编码。上层保证不会发生；
               真发生时跳过该字节会破坏往返，所以直接用最后一个符号顶替并不可取，
               这里选择保持前驱不变、跳过本字节，由上层的往返校验暴露问题。 */
            continue;
        }
        int cyc = fcv2_impl::cycleOf(i, len, rev);
        size_t base[MODEL_COUNT];
        d->cm.slotBases(prev, prev2, cyc, rv, base);

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
        prev2 = prev;
        prev = sym;
    }
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
    memset(d->symbolOf, -1, sizeof(d->symbolOf));
    uint32_t freq[TREE_CAP];
    for (int s = 0; s < alpha; s++) {
        int b = (int)d->rc.DecodeByte();
        int hi = (int)d->rc.DecodeByte();
        int lo = (int)d->rc.DecodeByte();
        d->byteOf[s] = b;
        d->symbolOf[b] = s;
        freq[s] = (uint32_t)((hi << 8) | lo);
    }
    d->alphaSize = alpha;
    /* 用码流里带来的量化频率重建，保证与编码端是同一棵树。 */
    d->tree.build(freq, alpha);
    d->cm.init(alpha);
    return 0;
}

int32_t coder_fcv2::decode_record(uint8_t* dst, uint32_t len, bool rev)
{
    fcv2_impl* d = impl.get();
    if (d->alphaSize == 0 || dst == nullptr || len == 0) {
        return 0;
    }
    int prev = d->alphaSize;
    int prev2 = d->alphaSize;
    int rv = rev ? 1 : 0;

    for (uint32_t i = 0; i < len; i++) {
        int cyc = fcv2_impl::cycleOf(i, len, rev);
        size_t base[MODEL_COUNT];
        d->cm.slotBases(prev, prev2, cyc, rv, base);

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
        prev2 = prev;
        prev = cur;
    }
    return (int32_t)len;
}
