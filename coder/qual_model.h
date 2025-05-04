#include <vector>
#include <stdint.h>
#include <memory.h>
#include "clr.h"
#include "../manager.h"

#define ABS(a)   ((a)>0?(a):-(a))
#define MIN(a, b) ((a)<(b)?(a):(b))
#define MAX(a, b) ((a)>(b)?(a):(b))

// Shrinking this to 1<<10 gives 2-3% smaller qualities, but 50% longer
#define QUAL_MAX_FREQ (1 << 16) - 32

typedef std::vector<std::pair<uint16_t, uint16_t>> SYMFREQ;

struct QUAL_MODEL
{
    enum
    {
        STEP = 8//18//16//12 , 18更好，20好像是效果最好的，24会变差 ,   TODO :根据数据类型，以及同一个数据动态设置STEP
    };

    QUAL_MODEL(){};

    virtual ~QUAL_MODEL() {}

    inline void reset(const SYMFREQ &sf);

    inline void init(const SYMFREQ &sf);

    inline void encodeSymbol(RangeCoder *rc, uint16_t sym);

    inline int encodeNearSymbol(RangeCoder *rc, uint16_t sym, int dist);

    inline uint16_t decodeSymbol(RangeCoder *rc);

    inline void free_freq();

    //protected:
    inline void normalize();

    uint32_t TotFreq; // Total frequency
    uint32_t BubCnt;  // Periodic counter for bubble sort step

    // Array of Symbols approximately sorted by Freq.
    struct SymFreqs
    {
        uint16_t Symbol;
        uint16_t Freq;
    } sentinel; // , F[nsym + 1]

    struct SymFreqs *F;
};

inline void QUAL_MODEL::init(const SYMFREQ &sf)
{
    safe_alloc(sf.size() + 1, struct SymFreqs, F);
    reset(sf);
}

inline void QUAL_MODEL::free_freq()
{
    free(F);
}

inline void QUAL_MODEL::reset(const SYMFREQ &sf)
{
    TotFreq = 0;
    for (size_t i = 0; i < sf.size(); i++)
    {
        F[i].Symbol = sf[i].first;
        F[i].Freq = sf[i].second;
        TotFreq += sf[i].second;
    }
    sentinel.Symbol = 0;
    sentinel.Freq = QUAL_MAX_FREQ; // Always first; simplifies sorting.
    BubCnt = 0;

    F[sf.size()].Freq = 0; // terminates normalize() loop. See below.
}

inline void QUAL_MODEL::normalize()
{
    /* Faster than F[i].Freq for 0 <= i < nsym */
    TotFreq = 0;
    for (SymFreqs *s = F; s->Freq; s++)
    {
        s->Freq -= s->Freq >> 1;
        TotFreq += s->Freq;
    }
}

inline void QUAL_MODEL::encodeSymbol(RangeCoder *rc, uint16_t sym)
{
    SymFreqs *s = F;
    uint32_t AccFreq = 0;

    while (s->Symbol != sym)
        AccFreq += s++->Freq;

    rc->Encode(AccFreq, s->Freq, TotFreq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > QUAL_MAX_FREQ)
        normalize();

    /* Keep approx sorted */ 
    // If the model is initialized with custom symbols (not in order like 1,2,3...) and frequencies, we can't sort, or decompression will fail. In practice, sorting has minimal impact on compression speed.
    // if (((++BubCnt & 15) == 0) && s[0].Freq > s[-1].Freq)
    // {
    //     SymFreqs t = s[0];
    //     s[0] = s[-1];
    //     s[-1] = t;
    // }
}

inline uint16_t QUAL_MODEL::decodeSymbol(RangeCoder *rc)
{
    SymFreqs *s = F;
    uint32_t freq = rc->GetFreq(TotFreq);
    uint32_t AccFreq;

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++)
        ;
    AccFreq -= s->Freq;

    rc->Decode(AccFreq, s->Freq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > QUAL_MAX_FREQ)
        normalize();

    /* Keep approx sorted */ 
    // If the model is initialized with custom symbols (not in order like 1,2,3...) and frequencies, we can't sort, or decompression will fail. In practice, sorting has minimal impact on compression speed.
    // if (((++BubCnt & 15) == 0) && s[0].Freq > s[-1].Freq)
    // {
    //     SymFreqs t = s[0];
    //     s[0] = s[-1];
    //     s[-1] = t;
    //     return t.Symbol;
    // }

    return s->Symbol;
}

struct QUAL_MODEL_ENGINE
{
public:
    QUAL_MODEL_ENGINE(const SYMFREQ &sf, const int64_t model_cnt)
    {
        int64_t i;
        sym_freq_initial.assign(sf.begin(), sf.end());
        ctx_cnt = model_cnt;
        
        // Replace problematic allocation with direct malloc
        model = static_cast<QUAL_MODEL*>(calloc(model_cnt, sizeof(QUAL_MODEL)));
        if (!model) {
            fprintf(stderr, "Error: Insufficient memory: need %" PRIu64 " MB\n",
                   (static_cast<uint64_t>(sizeof(QUAL_MODEL)) * model_cnt) >> 20);
            manage::instance().exit(ERR_MEM_NOENOUGH);
        }

        for (i = 0; i < model_cnt; i++)
            model[i].init(sf);
    }

    virtual ~QUAL_MODEL_ENGINE()
    {
        int64_t i;
        for (i = 0; i < ctx_cnt; i++)
            model[i].free_freq();
        free(model);
    }

    inline void encodeSymbol(RangeCoder *rc, uint16_t sym, const int64_t context)
    {
        return model[context].encodeSymbol(rc, sym);
    }

    inline uint16_t decodeSymbol(RangeCoder *rc, const int64_t context)
    {
        return model[context].decodeSymbol(rc);
    }

private:
    int64_t ctx_cnt;
    QUAL_MODEL *model;
    SYMFREQ sym_freq_initial; // Initialize which symbols the model has and their initial frequency values
    uint64_t tot_qual = 0; // Sum of quality values in current line
};