/*
 *--------------------------------------------------------------------------
 * A simple frequency model.
 *
 * This keeps a list of symbols and their frequencies, approximately
 * sorted by symbol frequency. We allow for a single symbol to periodically
 * move up the list when emitted, effectively doing a single step of
 * bubble sort periodically. This means it's largely the same complexity
 * irrespective of alphabet size.
 * It's more efficient on strongly biased distributions than random data.
 *
 * There is no escape symbol, so the model is tailored to relatively
 * stationary samples (although we do have occasional normalisation to
 * avoid frequency counters getting too high).
 *--------------------------------------------------------------------------
 */

// Shrinking this to 1<<10 gives 2-3% smaller qualities, but 50% longer

#ifndef _SIMPLE_MODEL_H_
#define _SIMPLE_MODEL_H_

#include "clr.h"

#define MAX_FREQ ((1<<16)-32)

template<int NSYM>
struct SIMPLE_MODEL {
    enum {
        STEP = 8
    };

    SIMPLE_MODEL();
    
    inline void reset(); 

    inline void encodeSymbol(RangeCoder *rc, uint16_t sym);
    inline void encodeSymbolOrder(RangeCoder *rc, uint16_t sym);
    inline void encodeSymbolNoUpdate(RangeCoder *rc, uint16_t sym);

    inline void updateModelAccFrecs(bool decode);

    inline uint16_t decodeSymbol(RangeCoder *rc);
    inline uint16_t decodeSymbolOrder(RangeCoder *rc);
    inline uint16_t decodeSymbolNoUpdate(RangeCoder *rc);

    void normalize();

    uint32_t TotFreq;  // Total frequency
    uint32_t BubCnt;   // Periodic counter for bubble sort step

    // Array of Symbols approximately sorted by Freq. 
    struct SymFreqs {
        uint16_t Symbol;
        uint16_t Freq;
        uint32_t AccFreq;
    } sentinel, F[NSYM + 1];
};


template<int NSYM>
SIMPLE_MODEL<NSYM>::SIMPLE_MODEL() {
    reset();
}

template<int NSYM>
void SIMPLE_MODEL<NSYM>::reset() {
    for (int i = 0; i < NSYM; i++) {
        F[i].Symbol = i;
        F[i].Freq = 1;
    }

    TotFreq = NSYM;
    sentinel.Symbol = 0;
    sentinel.Freq = MAX_FREQ; // Always first; simplifies sorting.
    BubCnt = 0;

    F[NSYM].Freq = 0; // terminates normalize() loop. See below.
}


template<int NSYM>
void SIMPLE_MODEL<NSYM>::normalize() {
    /* Faster than F[i].Freq for 0 <= i < NSYM */
    TotFreq = 0;
    for (SymFreqs *s = F; s->Freq; s++) {
        s->Freq -= s->Freq >> 1;
        TotFreq += s->Freq;
    }
}

template<int NSYM>
inline void SIMPLE_MODEL<NSYM>::encodeSymbol(RangeCoder *rc, uint16_t sym) {
    SymFreqs *s = F;
    uint32_t AccFreq = 0;

    while (s->Symbol != sym)
        AccFreq += s++->Freq;

    rc->Encode(AccFreq, s->Freq, TotFreq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > MAX_FREQ)
        normalize();
}

template<int NSYM>
inline void SIMPLE_MODEL<NSYM>::encodeSymbolOrder(RangeCoder *rc, uint16_t sym) {
    SymFreqs *s = F;
    uint32_t AccFreq = 0;
    // fprintf(stderr, "1=>s->Symbol: %u, sym: %u\n",s->Symbol, sym);
    while (s->Symbol != sym) {
        // fprintf(stderr, "2=>s->Symbol: %u, sym: %u\n",s->Symbol, sym);
        AccFreq += s++->Freq;
    }

    rc->Encode(AccFreq, s->Freq, TotFreq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > MAX_FREQ)
        normalize();

    /* Keep approx sorted */
    if (((++BubCnt & 15) == 0) && s[0].Freq > s[-1].Freq) {
        SymFreqs t = s[0];
        s[0] = s[-1];
        s[-1] = t;
    }
}

template<int NSYM>
inline void SIMPLE_MODEL<NSYM>::encodeSymbolNoUpdate(RangeCoder *rc, uint16_t sym) {
#ifdef __ORDER_SYMBOLS__
    SymFreqs *s = F;

    while (s->Symbol != sym)
        s++;

    rc->Encode(s->AccFreq, s->Freq, TotFreq);
#else
    rc->Encode(F[sym].AccFreq, F[sym].Freq, TotFreq);
#endif
}

template<int NSYM>
inline void SIMPLE_MODEL<NSYM>::updateModelAccFrecs(bool decode) {
#ifdef  __ORDER_SYMBOLS__
    qsort(F, NSYM, sizeof(SymFreqs), compare);
#endif
    if (decode) {
        F[0].AccFreq = F[0].Freq;
        for (int i = 1; i < NSYM; i++) {
            F[i].AccFreq = F[i - 1].AccFreq + F[i].Freq;
        }
    } else {
        F[0].AccFreq = 0;
        for (int i = 1; i < NSYM; i++) {
            F[i].AccFreq = F[i - 1].AccFreq + F[i - 1].Freq;
        }
    }
}

template<int NSYM>
inline uint16_t SIMPLE_MODEL<NSYM>::decodeSymbol(RangeCoder *rc) {
    SymFreqs *s = F;
    uint32_t freq = rc->GetFreq(TotFreq);
    uint32_t AccFreq;

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++);
    AccFreq -= s->Freq;

    rc->Decode(AccFreq, s->Freq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > MAX_FREQ)
        normalize();

    return s->Symbol;
}

template<int NSYM>
inline uint16_t SIMPLE_MODEL<NSYM>::decodeSymbolOrder(RangeCoder *rc) {
    SymFreqs *s = F;
    uint32_t freq = rc->GetFreq(TotFreq);
    uint32_t AccFreq;

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++);
    AccFreq -= s->Freq;

    rc->Decode(AccFreq, s->Freq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > MAX_FREQ)
        normalize();

    /* Keep approx sorted */
    if (((++BubCnt & 15) == 0) && s[0].Freq > s[-1].Freq) {
        SymFreqs t = s[0];
        s[0] = s[-1];
        s[-1] = t;
        return t.Symbol;
    }

    return s->Symbol;
}

template<int NSYM>
inline uint16_t SIMPLE_MODEL<NSYM>::decodeSymbolNoUpdate(RangeCoder *rc) {


    uint32_t freq = rc->GetFreq(TotFreq);
    SymFreqs *s = F;

    for (; s->AccFreq <= freq; s++);

    rc->Decode(s->AccFreq - s->Freq, s->Freq);

    return s->Symbol;
}

#endif