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

#include <limits>

#include "coder.h"
#include "fc/rangecoder.h"

namespace {

/* 哈夫曼树的节点数上限，同时也限制了字母表大小。 */
const int TREE_CAP = 64;

/* 权重更新的步长移位量，等效学习率 1/4096。 */
const int WEIGHT_LR_SHIFT = 12;

/* 参与混合的模型个数。m0..m5 为既有六档，m6 为策略 4 新增的 read 平均质量分档上下文。 */
const int MODEL_COUNT = 7;

/* m6 的 read 平均质量分档数（策略 4）。 */
const int QA_BINS = 4;

/* 上下文参数档位的合法范围，解码端读回码流头部时按同一规则归一，防止损坏的码流
   分配出过大的模型数组。各数组大小随 cfg 增长，必须给上界。 */
const int CFG_CYCLE_MAX_MIN = 32;
const int CFG_CYCLE_MAX_MAX = 128;
const int CFG_CYCLE_BUCKET_MAX = 32;
const int CFG_DELTA_MAX_MAX = 256;
const int CFG_DELTA_BUCKET_MAX = 16;
const int CFG_PREV_SHIFT_MAX = 3;

/*
 * 把外部传入的 cfg 归一成编码器实际使用的档位。两边（编码端、解码端）都必须调用，
 * 任何一条与这里不一致，编解码就会在上下文索引上错位。
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
        /* 关闭跳变上下文时把分档归一为 1，模型退化为"位置 + 树节点"，数组布局不变。 */
        cfg.deltaBucket = 1;
    }
    return cfg;
}

/*
 * 计算一条 read 的平均质量档（策略 4，fqzcomp 的 do_qa）。均值与存储顺序无关
 * （求和对整条质量串），按 phred 均值分 4 档。只有编码端用它把档位写进码流，
 * 解码端直接读回，所以这里的阈值只影响压缩率、不影响两端一致性。
 */
inline int qaBinOf(const uint8_t* qual, uint32_t len)
{
    uint64_t sum = 0;
    for (uint32_t i = 0; i < len; i++) {
        sum += qual[i];
    }
    int avg = (int)(sum / len);   /* 原始字节值，约 33..74（phred 0..41） */
    if (avg < 63) return 0;       /* phred < 30 */
    if (avg < 68) return 1;       /* 30..34 */
    if (avg < 71) return 2;       /* 35..37 */
    return 3;                     /* >= 38 */
}

/* 碱基上下文的取值个数：ACGTN + 未知（seq 为空或非碱基字符）。 */
const int BASE_STATES = 6;

/* 原始字节 -> 碱基符号（ACGTN -> 0..4），其余全部归入 5（未知）。initTables 里填充。 */
uint8_t g_baseSym[256];

/*
 * 模型快照与码流刻意分离：码流格式不能因先验模型功能而改变。下面的固定头让快照能在
 * 读入前判定是否仍对应当前的数组布局；只要任一尺寸常量变化，旧快照中的节点下标就不再
 * 有定义，必须拒绝，而不是以看似可用但实际错误的概率继续编码。
 */
const uint8_t MODEL_MAGIC[] = { 'f', 'c', 'v', '2', 'p', 'r', 'i', 'o', 'r' };
const uint16_t MODEL_FORMAT_VERSION = 7;  /* v7: read 平均质量档上下文 m6（useQa + QA_BINS） */

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
 * 字节读取器只在输入还剩足够空间时推进位置，因而截断快照不会越界读取，也不会让调用者
 * 看到只恢复了一半的模型。所有数字都逐字节按小端解析，避免宿主端序和结构体填充参与
 * 快照格式。
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
 * 七个粒度递增的上下文模型。
 *
 *   m0  只看树节点，最粗，任何时候样本都充足，起保底作用
 *   m1  位置 + 树节点
 *   m2  前一个符号 + 位置 + 树节点
 *   m3  前两个符号 + 位置档 + 链方向 + 树节点
 *   m4  当前碱基 + 位置 + 树节点
 *   m5  read 内质量跳变次数档 + 位置 + 树节点（策略 1）
 *   m6  read 平均质量分档 + 位置 + 树节点（策略 4）
 *
 * m4 用的碱基是"产生本质量值的那次测序循环"的碱基：正向链是位置 i，反向链按存储序
 * 取 len-1-i。质量值与碱基强相关（错配位、低复杂度区的质量系统性偏低），这维上下文
 * 与其他维度（前驱符号、链方向）基本正交，与 m3 不重叠。
 *
 * m5 的跳变次数是质量值随 read 单调变化（测序循环越靠后质量越低）的另一种刻画：
 * "到当前位置为止质量值一共变了多少次"反映这条 read 的质量是否干净，与局部前驱符号
 * （m2/m3）和位置（m1）都不重叠，且只在 read 内部累积、逐条重置。fqzcomp 用同一维
 * 度（state->delta）参与上下文，这里做成独立模型交给混合器，由权重决定它值多少。
 *
 * m6 是整条 read 的平均质量档（策略 4，fqzcomp 的 do_qa）：不同 read 的平均质量差异
 * 很大（实测 Q10~Q40 长尾），整条 read 处于"高质量档"还是"低质量档"是符号分布的
 * 强先验。档位由编码端算好写进码流（每记录 2 bit），解码端直接读回，上下文两侧一致。
 *
 * 各数组的尺寸由 Fcv2Cfg 决定，见 init()。
 */
struct ContextModel {
    int alphaSize;
    int prevStates;          /* 符号取值数加一，多出来的那个表示"没有前驱" */
    int prevQStates;         /* m3 用的量化前驱取值数（>> prevShift） */
    int cycleMax;
    int cycleBucket;
    int deltaMax;
    int deltaBucket;
    int prevShift;
    std::vector<Counter> m0, m1, m2, m3, m4, m5, m6;
    std::vector<int> weight; /* 混合权重，按 位置档 × 树节点 × 模型 组织 */

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
         * m3 用全符号取值太稀疏（本数据 39×39×16×2×64 ≈ 3.1M 槽，每块仅 ~9 次访问），
         * 学不出 order-2 结构。实测把前驱质量值 >>1 量化到 22 档（质量值≈符号 id），
         * 保留 55% 的 order-2 条件熵收益（3.66 vs 3.53 bits），槽数却降到 1/3。
         * 档位随 prevShift 可调：字母表大时继续右移更粗，字母表小时可以不量化。
         */
        prevQStates = (alpha >> prevShift) + 2;
        m0.assign((size_t)TREE_CAP, COUNTER_INIT);
        m1.assign((size_t)cycleMax * TREE_CAP, COUNTER_INIT);
        m2.assign((size_t)prevStates * cycleMax * TREE_CAP, COUNTER_INIT);
        m3.assign((size_t)prevQStates * prevQStates * cycleBucket * 2 * TREE_CAP, COUNTER_INIT);
        m4.assign((size_t)BASE_STATES * cycleMax * TREE_CAP, COUNTER_INIT);
        m5.assign((size_t)deltaBucket * cycleMax * TREE_CAP, COUNTER_INIT);
        m6.assign((size_t)QA_BINS * cycleMax * TREE_CAP, COUNTER_INIT);
        /* 权重初值取 1<<14，即定点 0.25，七个模型合起来接近简单平均。 */
        weight.assign((size_t)cycleBucket * TREE_CAP * MODEL_COUNT, 1 << 14);
    }

    inline void slotBases(int prev, int prev2, int cyc, int rev, int baseSym,
                          int delta, int qa, size_t* base) const
    {
        int bucket = cyc * cycleBucket / cycleMax;
        if (bucket >= cycleBucket) bucket = cycleBucket - 1;
        int qp = (prev >> prevShift);      /* 量化前驱质量值 */
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
 * 实现体。放在匿名命名空间之外是因为头文件里前置声明了它。
 */
class fcv2_impl {
public:
    coder_io* io;
    ContextModel cm;
    HuffTree tree;
    RangeCoder rc;

    /*
     * 链方向自己的概率模型。
     *
     * 链方向随每条记录写进码流，而不是由调用方在解压时再提供一次。这样编码器完全
     * 自包含：解压侧不需要为了解 QUAL 去关心 FLAG 字段解到哪一步了，也不必维护一份
     * 逐记录的链方向数组。
     *
     * 代价是每条记录约一个比特。实测反向链占 49.9%，接近最大熵，压不动，一百万条
     * 记录约 125 KB，对 90 MB 的质量值是 0.139%。而链方向带来的收益是 0.372 个
     * 百分点，扣掉这部分仍然净赚 0.233 个百分点，值得。
     *
     * 仍然用一个自适应计数器而不是固定的 0.5，是因为按坐标排序的 SAM 里正反链未必
     * 严格各半，模型能自己适应实际比例。
     */
    Counter revCounter;

    /* 相邻重复 read 去重（策略 3）：dup 标记自己的自适应概率，随记录逐条维护。 */
    Counter dupCounter;
    /* 上一条记录的质量串（存储序），用于比对与解码端还原；跨记录保留。 */
    std::vector<uint8_t> prevQual;
    uint32_t prevQualLen = 0;

    Fcv2Cfg cfg;             /* 归一后的上下文参数档位 */

    int alphaSize;
    int symbolOf[256];      /* 原始字节 -> 内部符号编号 */
    int byteOf[TREE_CAP];   /* 内部符号编号 -> 原始字节 */
    uint16_t quantFreq[TREE_CAP]; /* 量化后的频率，随码流写出供解码端重建同一棵树 */
    bool encodeStarted;
    bool flushed;
    bool modelLoaded;

    explicit fcv2_impl(coder_io* ioPtr, const Fcv2Cfg& cfgIn)
        : io(ioPtr), cfg(normalizeCfg(cfgIn)), revCounter(COUNTER_INIT),
          dupCounter(COUNTER_INIT), alphaSize(0),
          encodeStarted(false), flushed(false), modelLoaded(false)
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
         * 先验是跨块训练产物，其计数器数组按训练时的档位布局。加载方（可能是压缩端或
         * 解码端）采用先验自身携带的档位作为有效档位，而不是要求先验与构造时的档位
         * 一致：先验总是与它一起使用的那条码流同档位（训练与压缩用同一组参数），
         * 这样任何一端都能安全地"先验决定了档位"。若与码流头部不一致，begin_decode
         * 会以档位不符拒绝，见其 modelLoaded 分支。
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
        prevQualLen = 0;   /* 先验只带模型状态，不带上一条质量串 */
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
         * 上下文参数档位写进码流头部。模型数组的布局由它决定，解码端必须读回同一组
         * 参数才能对上上下文索引；这也让每个数据块可以独立携带自己的档位选择。
         * 见 normalizeCfg()：两端按同一规则归一，保证任何合法取值两侧一致。
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
     * 相邻重复 read 去重（策略 3，fqzcomp 的 do_dedup）：与上一条记录逐字节比对，
     * 完全相同就写 1 bit 并把整条质量串跳过。命中重复时本记录的上下文状态不更新
     * （prev/prev2/delta 都在记录开头重置），下一条记录照常从零开始，两端一致。
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
     * read 平均质量档（策略 4）：档位写进码流，作为 m6 上下文。固定 2 bit/记录，
     * 解码端读回同一档位。
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
    int delta = 0;   /* read 内质量跳变次数，随记录重置 */

    for (uint32_t i = 0; i < len; i++) {
        int sym = d->symbolOf[qual[i]];
        if (sym < 0) {
            /* 出现了统计阶段没见过的质量值，无法编码。上层保证不会发生；
               真发生时跳过该字节会破坏往返，所以直接用最后一个符号顶替并不可取，
               这里选择保持前驱不变、跳过本字节，由上层的往返校验暴露问题。 */
            continue;
        }
        /*
         * 碱基上下文取"产生本质量值的循环"对应的碱基：正向链 i、反向链 len-1-i，
         * 与 cycleOf 用同一套映射，保证两端一致。seq 长度不足（异常数据）时该位置
         * 归入"未知"档，绝不对 seq 越界读取。
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
         * 跳变次数在本符号编码后才更新，因此 m5 的上下文统计的是"当前符号之前"的变化
         * 次数，与 fqzcomp 的 delta 用法一致。第一个符号没有前驱，不算跳变。
         */
        if (prev != d->alphaSize && sym != prev) {
            delta++;
        }
        prev2 = prev;
        prev = sym;
    }

    /* 非重复记录：记录本条质量串，供下一条做去重比对。 */
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
    /* 读回编码端写入的上下文参数档位，与 encode_record 的写出顺序严格一致。 */
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
        return -1;   /* 损坏码流：档位超出合法范围，禁止据此分配模型 */
    }
    streamCfg = normalizeCfg(streamCfg);
    /*
     * 先验来自跨块训练，其模型布局由训练时的档位决定；如果与码流头部不一致，计数器
     * 数组对不上，只能拒绝（压缩侧的同款不一致会让 loadModel 返回 false）。
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
     * 先验中的计数器按树节点编号索引，不能拿它去解另一棵树的码流。即使两个字母表含有
     * 同一批字节，只要频率改变导致合并顺序变化，节点编号也可能已经代表另一个判断；因此
     * 载入先验后必须逐项核对码流头，并保留已载入的学习状态，而非再次 init 覆盖它。
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
    /* 用码流里带来的量化频率重建，保证与编码端是同一棵树。 */
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

    /* 链方向由码流自带，见 fcv2_impl::revCounter 的说明。 */
    int p0 = counterProb(d->revCounter);
    if (p0 < 1)    p0 = 1;
    if (p0 > 4095) p0 = 4095;
    int rv = d->rc.DecodeBit<12>(p0);
    counterUpdate(d->revCounter, !rv);
    bool rev = (rv != 0);

    /* 相邻重复 read 去重：dup=1 时直接拷贝上一条解码出的质量串，跳过本记录。 */
    if (d->cfg.useDedup) {
        int pDup = counterProb(d->dupCounter);
        if (pDup < 1)    pDup = 1;
        if (pDup > 4095) pDup = 4095;
        int isDup = d->rc.DecodeBit<12>(pDup);
        counterUpdate(d->dupCounter, !isDup);
        if (isDup) {
            if (d->prevQualLen < len) {
                /* 上一条质量串缺失或过短，属于损坏流；用 '!' 兜底避免越界读。 */
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
    int delta = 0;   /* 与编码侧同序维护的 read 内质量跳变次数 */

    /* read 平均质量档：读回编码端写下的 2 bit，作为 m6 上下文。 */
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
