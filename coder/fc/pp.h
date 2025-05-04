#ifndef _PP_H_
#define _PP_H_

#include "fc_header.h"
#include "../../manager.h"

#define FC_MATCH 	0xf2

#if defined (FC_UNALIGN) && defined (__x86_64__)

template<class T> int fc_pp_2(const unsigned char * RESTRICT input, const unsigned char * inputEnd, unsigned char * RESTRICT output, unsigned char * outputEnd, int * RESTRICT lookup, int mask)
{
    const unsigned char *   inputStart      = input;
    const unsigned char *   inputMinLenEnd  = inputEnd - sizeof(T) - 32;

    const unsigned char *   outputStart     = output;
    const unsigned char *   outputEOB       = outputEnd - 8;

    for (int i = 0; i < 4; ++i) { *output++ = *input++; }

    {
        while ((input < inputMinLenEnd) && (output < outputEOB))
        {
            unsigned long long next8 = *(unsigned long long *)(input - 4); *(unsigned int *)(output) = (unsigned int)(next8 >> 32); next8 = fc_byteswap_uint64(next8);

            int value;
            {
                const unsigned int index0 = (((next8 >> (4 * 8)) >> 15) ^ (next8 >> (4 * 8)) ^ ((next8 >> (4 * 8)) >> 3)) & mask; value = lookup[index0]; lookup[index0] = (int)(input - inputStart + 0); 
                if (value > 0 && (*(T *)(input + 0) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND1;
                if (value > 0 && ((unsigned char)(next8 >> 3 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND1;

                const unsigned int index1 = (((next8 >> (3 * 8)) >> 15) ^ (next8 >> (3 * 8)) ^ ((next8 >> (3 * 8)) >> 3)) & mask; value = lookup[index1]; lookup[index1] = (int)(input - inputStart + 1); 
                if (value > 0 && (*(T *)(input + 1) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND2;
                if (value > 0 && ((unsigned char)(next8 >> 2 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND2;

                const unsigned int index2 = (((next8 >> (2 * 8)) >> 15) ^ (next8 >> (2 * 8)) ^ ((next8 >> (2 * 8)) >> 3)) & mask; value = lookup[index2]; lookup[index2] = (int)(input - inputStart + 2); 
                if (value > 0 && (*(T *)(input + 2) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND3;
                if (value > 0 && ((unsigned char)(next8 >> 1 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND3;

                const unsigned int index3 = (((next8 >> (1 * 8)) >> 15) ^ (next8 >> (1 * 8)) ^ ((next8 >> (1 * 8)) >> 3)) & mask; value = lookup[index3]; lookup[index3] = (int)(input - inputStart + 3); 
                if (value > 0 && (*(T *)(input + 3) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND4;
                if (value > 0 && ((unsigned char)(next8 >> 0 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND4;

                input += 4; output += 4;

                continue;
            }

FC_PP_GOOD_MATCH_FOUND4:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND3:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND2:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND1:

            {
                const unsigned char * RESTRICT reference = inputStart + value;

                long long len = sizeof(T);

                for (; input + len < inputMinLenEnd; len += sizeof(unsigned long long))
                {
                    unsigned long long m;
                    if ((m = (*(unsigned long long *)(input + len)) ^ *(unsigned long long *)(reference + len)) != 0) 
                    {
                        len += fc_bit_scan_forward64(m) / 8; break;
                    }
                }

                input += len; len -= sizeof(T);

                *output++ = FC_MATCH; while (len >= 254) { len -= 254; *output++ = 254; if (output >= outputEOB) break; } *output++ = (unsigned char)(len); 
            
                continue;
            }

FC_PP_BAD_MATCH_FOUND4:
            input += 4; output += 4; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND3:
            input += 3; output += 3; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND2:
            input += 2; output += 2; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND1:
            input += 1; output += 1; *output++ = 255; continue;
        }
    }
    
    {
        unsigned int context = input[-1] | (input[-2] << 8) | (input[-3] << 16) | (input[-4] << 24);

        while ((input < inputEnd) && (output < outputEOB))
        {
            unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
            int value = lookup[index]; lookup[index] = (int)(input - inputStart);

            unsigned char next = *output++ = *input++; context = (context << 8) | next;
            if (next == FC_MATCH && value > 0) *output++ = 255;
        }
    }

    return (output >= outputEOB) ? FC_FAILED_ZIP : (int)(output - outputStart);
}

template<class T> int fc_pp_32x(const unsigned char * RESTRICT input, const unsigned char * inputEnd, unsigned char * RESTRICT output, unsigned char * outputEnd, int * RESTRICT lookup, int mask)
{
    const unsigned char *   inputStart      = input;
    const unsigned char *   inputMinLenEnd  = inputEnd - sizeof(T) - sizeof(T) - 32;

    const unsigned char *   outputStart     = output;
    const unsigned char *   outputEOB       = outputEnd - 8;

    for (int i = 0; i < 4; ++i) { *output++ = *input++; }

    {
        while ((input < inputMinLenEnd) && (output < outputEOB))
        {
            unsigned long long next8 = *(unsigned long long *)(input - 4); *(unsigned int *)(output) = (unsigned int)(next8 >> 32); next8 = fc_byteswap_uint64(next8);

            int value;
            {
                const unsigned int index0 = (((next8 >> (4 * 8)) >> 15) ^ (next8 >> (4 * 8)) ^ ((next8 >> (4 * 8)) >> 3)) & mask; value = lookup[index0]; lookup[index0] = (int)(input - inputStart + 0); 
                if (value > 0 && (*(T *)(input + sizeof(T) + 0) == *(T *)(inputStart + value + sizeof(T))) && (*(T *)(input + 0) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND1;
                if (value > 0 && ((unsigned char)(next8 >> 3 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND1;

                const unsigned int index1 = (((next8 >> (3 * 8)) >> 15) ^ (next8 >> (3 * 8)) ^ ((next8 >> (3 * 8)) >> 3)) & mask; value = lookup[index1]; lookup[index1] = (int)(input - inputStart + 1); 
                if (value > 0 && (*(T *)(input + sizeof(T) + 1) == *(T *)(inputStart + value + sizeof(T))) && (*(T *)(input + 1) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND2;
                if (value > 0 && ((unsigned char)(next8 >> 2 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND2;

                const unsigned int index2 = (((next8 >> (2 * 8)) >> 15) ^ (next8 >> (2 * 8)) ^ ((next8 >> (2 * 8)) >> 3)) & mask; value = lookup[index2]; lookup[index2] = (int)(input - inputStart + 2); 
                if (value > 0 && (*(T *)(input + sizeof(T) + 2) == *(T *)(inputStart + value + sizeof(T))) && (*(T *)(input + 2) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND3;
                if (value > 0 && ((unsigned char)(next8 >> 1 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND3;

                const unsigned int index3 = (((next8 >> (1 * 8)) >> 15) ^ (next8 >> (1 * 8)) ^ ((next8 >> (1 * 8)) >> 3)) & mask; value = lookup[index3]; lookup[index3] = (int)(input - inputStart + 3); 
                if (value > 0 && (*(T *)(input + sizeof(T) + 3) == *(T *)(inputStart + value + sizeof(T))) && (*(T *)(input + 3) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND4;
                if (value > 0 && ((unsigned char)(next8 >> 0 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND4;

                input += 4; output += 4;

                continue;
            }

FC_PP_GOOD_MATCH_FOUND4:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND3:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND2:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND1:

            {
                const unsigned char * RESTRICT reference = inputStart + value;

                long long len = sizeof(T) + sizeof(T);

                for (; input + len < inputMinLenEnd; len += sizeof(unsigned long long))
                {
                    unsigned long long m;
                    if ((m = (*(unsigned long long *)(input + len)) ^ *(unsigned long long *)(reference + len)) != 0) 
                    {
                        len += fc_bit_scan_forward64(m) / 8; break;
                    }
                }

                input += len; len -= sizeof(T) + sizeof(T);

                *output++ = FC_MATCH; while (len >= 254) { len -= 254; *output++ = 254; if (output >= outputEOB) break; } *output++ = (unsigned char)(len); 
            
                continue;
            }

FC_PP_BAD_MATCH_FOUND4:
            input += 4; output += 4; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND3:
            input += 3; output += 3; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND2:
            input += 2; output += 2; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND1:
            input += 1; output += 1; *output++ = 255; continue;
        }
    }
    
    {
        unsigned int context = input[-1] | (input[-2] << 8) | (input[-3] << 16) | (input[-4] << 24);

        while ((input < inputEnd) && (output < outputEOB))
        {
            unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
            int value = lookup[index]; lookup[index] = (int)(input - inputStart);

            unsigned char next = *output++ = *input++; context = (context << 8) | next;
            if (next == FC_MATCH && value > 0) *output++ = 255;
        }
    }

    return (output >= outputEOB) ? FC_FAILED_ZIP : (int)(output - outputStart);
}

template<class T> int fc_pp_4(const unsigned char * RESTRICT input, const unsigned char * inputEnd, unsigned char * RESTRICT output, unsigned char * outputEnd, int * RESTRICT lookup, int mask, int minLen)
{
    const unsigned char *   inputStart      = input;
    const unsigned char *   inputMinLenEnd  = inputEnd - sizeof(T) - sizeof(T) - 32;

    const unsigned char *   outputStart     = output;
    const unsigned char *   outputEOB       = outputEnd - 8;

    for (int i = 0; i < 4; ++i) { *output++ = *input++; }

    {
        while ((input < inputMinLenEnd) && (output < outputEOB))
        {
            unsigned long long next8 = *(unsigned long long *)(input - 4); *(unsigned int *)(output) = (unsigned int)(next8 >> 32); next8 = fc_byteswap_uint64(next8);

            int value;
            {
                const unsigned int index0 = (((next8 >> (4 * 8)) >> 15) ^ (next8 >> (4 * 8)) ^ ((next8 >> (4 * 8)) >> 3)) & mask; value = lookup[index0]; lookup[index0] = (int)(input - inputStart + 0); 
                if (value > 0 && (*(T *)(input + minLen - sizeof(T) + 0) == *(T *)(inputStart + value + minLen - sizeof(T))) && (*(T *)(input + 0) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND1;
                if (value > 0 && ((unsigned char)(next8 >> 3 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND1;

                const unsigned int index1 = (((next8 >> (3 * 8)) >> 15) ^ (next8 >> (3 * 8)) ^ ((next8 >> (3 * 8)) >> 3)) & mask; value = lookup[index1]; lookup[index1] = (int)(input - inputStart + 1); 
                if (value > 0 && (*(T *)(input + minLen - sizeof(T) + 1) == *(T *)(inputStart + value + minLen - sizeof(T))) && (*(T *)(input + 1) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND2;
                if (value > 0 && ((unsigned char)(next8 >> 2 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND2;

                const unsigned int index2 = (((next8 >> (2 * 8)) >> 15) ^ (next8 >> (2 * 8)) ^ ((next8 >> (2 * 8)) >> 3)) & mask; value = lookup[index2]; lookup[index2] = (int)(input - inputStart + 2); 
                if (value > 0 && (*(T *)(input + minLen - sizeof(T) + 2) == *(T *)(inputStart + value + minLen - sizeof(T))) && (*(T *)(input + 2) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND3;
                if (value > 0 && ((unsigned char)(next8 >> 1 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND3;

                const unsigned int index3 = (((next8 >> (1 * 8)) >> 15) ^ (next8 >> (1 * 8)) ^ ((next8 >> (1 * 8)) >> 3)) & mask; value = lookup[index3]; lookup[index3] = (int)(input - inputStart + 3); 
                if (value > 0 && (*(T *)(input + minLen - sizeof(T) + 3) == *(T *)(inputStart + value + minLen - sizeof(T))) && (*(T *)(input + 3) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND4;
                if (value > 0 && ((unsigned char)(next8 >> 0 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND4;

                input += 4; output += 4;

                continue;
            }

FC_PP_GOOD_MATCH_FOUND4:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND3:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND2:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND1:

            {
                const unsigned char * RESTRICT reference = inputStart + value;

                long long len = minLen;

                for (; input + len < inputMinLenEnd; len += sizeof(unsigned long long))
                {
                    unsigned long long m;
                    if ((m = (*(unsigned long long *)(input + len)) ^ *(unsigned long long *)(reference + len)) != 0) 
                    {
                        len += fc_bit_scan_forward64(m) / 8; break;
                    }
                }

                input += len; len -= minLen;

                *output++ = FC_MATCH; while (len >= 254) { len -= 254; *output++ = 254; if (output >= outputEOB) break; } *output++ = (unsigned char)(len); 
            
                continue;
            }

FC_PP_BAD_MATCH_FOUND4:
            input += 4; output += 4; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND3:
            input += 3; output += 3; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND2:
            input += 2; output += 2; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND1:
            input += 1; output += 1; *output++ = 255; continue;
        }
    }
    
    {
        unsigned int context = input[-1] | (input[-2] << 8) | (input[-3] << 16) | (input[-4] << 24);

        while ((input < inputEnd) && (output < outputEOB))
        {
            unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
            int value = lookup[index]; lookup[index] = (int)(input - inputStart);

            unsigned char next = *output++ = *input++; context = (context << 8) | next;
            if (next == FC_MATCH && value > 0) *output++ = 255;
        }
    }

    return (output >= outputEOB) ? FC_FAILED_ZIP : (int)(output - outputStart);
}

#define NEXT8_CONTEXT_ID(n) ((next8 >> (n * 8)) ^ ((next8 >> (n * 8)) >> 1) ^ ((next8 >> (n * 8)) >> 2) ^ ((next8 >> (n * 8)) >> 3))

template<class T> int fc_pp_1(const unsigned char * RESTRICT input, const unsigned char * inputEnd, unsigned char * RESTRICT output, unsigned char * outputEnd, int * RESTRICT lookup, int mask, int minLen)
{
    const unsigned char *   inputStart  = input;
    const unsigned char *   outputStart = output;
    const unsigned char *   outputEOB   = outputEnd - 8;

    const unsigned char * heuristic      = input;
    const unsigned char * inputMinLenEnd = inputEnd - minLen - 32;

    for (int i = 0; i < 4; ++i) { *output++ = *input++; }

    {
        while ((input < inputMinLenEnd) && (output < outputEOB))
        {
            unsigned long long next8 = *(unsigned long long *)(input - 4); *(unsigned int *)(output) = (unsigned int)(next8 >> 32); next8 = fc_byteswap_uint64(next8);

            int value;
            {
                // const unsigned int index0 = (((next8 >> (4 * 8)) >> r1shift) ^ ((next8 >> (4 * 8)) >> 5)  ^ ((next8 >> (4 * 8)) >> 2)^ (next8 >> (4 * 8)) ^ ((next8 >> (4 * 8)) >> r2shift)) & mask; value = lookup[index0]; lookup[index0] = (int)(input - inputStart + 0); 
                const unsigned int index0 = NEXT8_CONTEXT_ID(4) & mask; value = lookup[index0]; lookup[index0] = (int)(input - inputStart + 0); 
                if (value > 0 && input > heuristic && (*(T *)(input + minLen - sizeof(T) + 0) == *(T *)(inputStart + value + minLen - sizeof(T))) && (*(T *)(input + 0) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND1;
                if (value > 0 && ((unsigned char)(next8 >> 3 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND1;

                const unsigned int index1 = NEXT8_CONTEXT_ID(3) & mask; value = lookup[index1]; lookup[index1] = (int)(input - inputStart + 1); 
                if (value > 0 && input > heuristic && (*(T *)(input + minLen - sizeof(T) + 1) == *(T *)(inputStart + value + minLen - sizeof(T))) && (*(T *)(input + 1) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND2;
                if (value > 0 && ((unsigned char)(next8 >> 2 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND2;

                const unsigned int index2 = NEXT8_CONTEXT_ID(2) & mask; value = lookup[index2]; lookup[index2] = (int)(input - inputStart + 2); 
                if (value > 0 && input > heuristic && (*(T *)(input + minLen - sizeof(T) + 2) == *(T *)(inputStart + value + minLen - sizeof(T))) && (*(T *)(input + 2) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND3;
                if (value > 0 && ((unsigned char)(next8 >> 1 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND3;

                const unsigned int index3 = NEXT8_CONTEXT_ID(1) & mask; value = lookup[index3]; lookup[index3] = (int)(input - inputStart + 3); 
                if (value > 0 && input > heuristic && (*(T *)(input + minLen - sizeof(T) + 3) == *(T *)(inputStart + value + minLen - sizeof(T))) && (*(T *)(input + 3) == *(T *)(inputStart + value))) goto FC_PP_GOOD_MATCH_FOUND4;
                if (value > 0 && ((unsigned char)(next8 >> 0 * 8) == FC_MATCH)) goto FC_PP_BAD_MATCH_FOUND4;

                input += 4; output += 4;
            
                continue;
            }

FC_PP_GOOD_MATCH_FOUND4:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND3:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND2:
            input += 1; output += 1;
FC_PP_GOOD_MATCH_FOUND1:

            {
                const unsigned char * RESTRICT reference = inputStart + value;

                long long len = sizeof(T);

                for (; input + len < inputMinLenEnd; len += sizeof(unsigned long long))
                {
                    unsigned long long m;
                    if ((m = (*(unsigned long long *)(input + len)) ^ *(unsigned long long *)(reference + len)) != 0) 
                    {
                        len += fc_bit_scan_forward64(m) / 8; break;
                    }
                }

                if (len < minLen) { heuristic = input + len; goto FC_PP_MATCH_NOT_FOUND; }

                input += len; len -= minLen;

                *output++ = FC_MATCH; while (len >= 254) { len -= 254; *output++ = 254; if (output >= outputEOB) break; } *output++ = (unsigned char)(len); 
            
                continue;
            }

FC_PP_MATCH_NOT_FOUND:
            if ((*output++ = *input++) == FC_MATCH) { *output++ = 255; }

            continue;

FC_PP_BAD_MATCH_FOUND4:
            input += 4; output += 4; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND3:
            input += 3; output += 3; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND2:
            input += 2; output += 2; *output++ = 255; continue;
FC_PP_BAD_MATCH_FOUND1:
            input += 1; output += 1; *output++ = 255; continue;
        }        
    }
    
    {
        unsigned int context = input[-1] | (input[-2] << 8) | (input[-3] << 16) | (input[-4] << 24);

        while ((input < inputEnd) && (output < outputEOB))
        {
            unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
            int value = lookup[index]; lookup[index] = (int)(input - inputStart);

            unsigned char next = *output++ = *input++; context = (context << 8) | next;
            if (next == FC_MATCH && value > 0) *output++ = 255;
        }
    }

    return (output >= outputEOB) ? FC_FAILED_ZIP : (int)(output - outputStart);
}

#endif

static int fc_pp_encode_generic(const unsigned char * RESTRICT input, const unsigned char * inputEnd, unsigned char * RESTRICT output, unsigned char * outputEnd, int * RESTRICT lookup, int mask, int minLen)
{
    const unsigned char *   inputStart  = input;
    const unsigned char *   outputStart = output;
    const unsigned char *   outputEOB   = outputEnd - 8;

    const unsigned char * heuristic      = input;
    const unsigned char * inputMinLenEnd = inputEnd - minLen - 32;

    for (int i = 0; i < 4; ++i) { *output++ = *input++; }

    {
        unsigned int context = input[-1] | (input[-2] << 8) | (input[-3] << 16) | (input[-4] << 24);

        while ((input < inputMinLenEnd) && (output < outputEOB))
        {
            unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
            int value = lookup[index]; lookup[index] = (int)(input - inputStart);
            if (value > 0)
            {
                const unsigned char * RESTRICT reference = inputStart + value;
#if defined(FC_UNALIGN)
                if ((*(unsigned int *)(input + minLen - 4) == *(unsigned int *)(reference + minLen - 4)) && (*(unsigned int *)(input) == *(unsigned int *)(reference)))
#else
                if ((memcmp(input + minLen - 4, reference + minLen - 4, sizeof(unsigned int)) == 0) && (memcmp(input, reference, sizeof(unsigned int)) == 0))
#endif
                {
                    if ((heuristic > input) && (*(unsigned int *)heuristic != *(unsigned int *)(reference + (heuristic - input))))
                    {
                        goto FC_PP_MATCH_NOT_FOUND;
                    }

                    int len = 4;
                    for (; input + len < inputMinLenEnd; len += sizeof(unsigned int))
                    {
                        if (*(unsigned int *)(input + len) != *(unsigned int *)(reference + len)) break;
                    }

                    if (len < minLen)
                    {
                        if (heuristic < input + len) heuristic = input + len;
                        goto FC_PP_MATCH_NOT_FOUND;
                    }

#if defined(FC_UNALIGN)
                    len += sizeof(unsigned short) * (*(unsigned short *)(input + len) == *(unsigned short *)(reference + len));
                    len += sizeof(unsigned char ) * (*(unsigned char  *)(input + len) == *(unsigned char  *)(reference + len));
#else
                    len += input[len] == reference[len];
                    len += input[len] == reference[len];
                    len += input[len] == reference[len];
#endif

                    input += len; context = input[-1] | (input[-2] << 8) | (input[-3] << 16) | (input[-4] << 24);

                    *output++ = FC_MATCH;

                    len -= minLen; while (len >= 254) { len -= 254; *output++ = 254; if (output >= outputEOB) break; }

                    *output++ = (unsigned char)(len);
                }
                else
                {

FC_PP_MATCH_NOT_FOUND:
                    unsigned char next = *output++ = *input++; context = (context << 8) | next;
                    if (next == FC_MATCH) *output++ = 255;
                }
            }
            else
            {
                context = (context << 8) | (*output++ = *input++);
            }
        }
    }
    
    {
        unsigned int context = input[-1] | (input[-2] << 8) | (input[-3] << 16) | (input[-4] << 24);

        while ((input < inputEnd) && (output < outputEOB))
        {
            unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
            int value = lookup[index]; lookup[index] = (int)(input - inputStart);

            unsigned char next = *output++ = *input++; context = (context << 8) | next;
            if (next == FC_MATCH && value > 0) *output++ = 255;
        }
    }

    return (output >= outputEOB) ? FC_FAILED_ZIP : (int)(output - outputStart);
}


static int fc_preprocess(const unsigned char * input, const unsigned char * inputEnd, unsigned char * output, unsigned char * outputEnd, int hashSize, int minLen)
{
    if (inputEnd - input - minLen < 32)
    {
        return FC_FAILED_ZIP;
    }

    int result = FC_LACK_OF_MEMORY;
    if (int * lookup = (int *)calloc(1, (int)(1 << hashSize) * sizeof(int))) 
    {
#if defined (FC_UNALIGN) && defined (__x86_64__)
        result = (minLen == 1 * (int)sizeof(unsigned int      ) && result == FC_LACK_OF_MEMORY) ? fc_pp_2  <unsigned int      >(input, inputEnd, output, outputEnd, lookup, (int)(1 << hashSize) - 1) : result;
        result = (minLen == 1 * (int)sizeof(unsigned long long) && result == FC_LACK_OF_MEMORY) ? fc_pp_2  <unsigned long long>(input, inputEnd, output, outputEnd, lookup, (int)(1 << hashSize) - 1) : result;
        result = (minLen == 2 * (int)sizeof(unsigned long long) && result == FC_LACK_OF_MEMORY) ? fc_pp_32x<unsigned long long>(input, inputEnd, output, outputEnd, lookup, (int)(1 << hashSize) - 1) : result;
        result = (minLen <= 2 * (int)sizeof(unsigned int      ) && result == FC_LACK_OF_MEMORY) ? fc_pp_4 <unsigned int      >(input, inputEnd, output, outputEnd, lookup, (int)(1 << hashSize) - 1, minLen) : result;
        result = (minLen <= 2 * (int)sizeof(unsigned long long) && result == FC_LACK_OF_MEMORY) ? fc_pp_4 <unsigned long long>(input, inputEnd, output, outputEnd, lookup, (int)(1 << hashSize) - 1, minLen) : result;
        
        result = result == FC_LACK_OF_MEMORY ? fc_pp_1<unsigned long long>(input, inputEnd, output, outputEnd, lookup, (int)(1 << hashSize) - 1, minLen) : result;
#endif

        result = result == FC_LACK_OF_MEMORY ? fc_pp_encode_generic(input, inputEnd, output, outputEnd, lookup, (int)(1 << hashSize) - 1, minLen) : result;

        free(lookup);
    }

    return result;
}


static int fc_unpreprocess(const unsigned char * RESTRICT input, const unsigned char * inputEnd, unsigned char * RESTRICT output, int hashSize, int minLen)
{
    if (inputEnd - input < 4)
    {
        return FC_INVALID_END;
    }

    if (int * RESTRICT lookup = (int *)calloc(1, ((int)(1 << hashSize) * sizeof(int))))
    {
        unsigned int            mask        = (int)(1 << hashSize) - 1;
        const unsigned char *   outputStart = output;

        for (int i = 0; i < 4; ++i) { *output++ = *input++; }

#if defined (FC_UNALIGN) && defined (__x86_64__)
        if (hashSize <= 17)
        {
            unsigned int prev4 = *(unsigned int *)(output - 4);

            while (input < inputEnd - 8)
            {
                unsigned int next4          = *(unsigned int *)(output) = *(unsigned int *)(input);
                unsigned long long next8    = fc_byteswap_uint64(((unsigned long long)next4 << 32) | prev4);

                int value;
                {
                    // const unsigned int index0 = (((next8 >> (4 * 8)) >> r1shift) ^ (next8 >> (4 * 8)) ^ ((next8 >> (4 * 8)) >> r2shift)) & mask;
                    const unsigned int index0 = NEXT8_CONTEXT_ID(4) & mask;
                    value = lookup[index0]; lookup[index0] = (int)(output - outputStart + 0); if (((unsigned char)(next8 >> 3 * 8) == FC_MATCH) && (value > 0)) goto FC_PP_MATCH_FOUND1;

                    const unsigned int index1 = NEXT8_CONTEXT_ID(3) & mask;
                    value = lookup[index1]; lookup[index1] = (int)(output - outputStart + 1); if (((unsigned char)(next8 >> 2 * 8) == FC_MATCH) && (value > 0)) goto FC_PP_MATCH_FOUND2;

                    const unsigned int index2 = NEXT8_CONTEXT_ID(2) & mask;
                    value = lookup[index2]; lookup[index2] = (int)(output - outputStart + 2); if (((unsigned char)(next8 >> 1 * 8) == FC_MATCH) && (value > 0)) goto FC_PP_MATCH_FOUND3;

                    const unsigned int index3 = NEXT8_CONTEXT_ID(1) & mask;
                    value = lookup[index3]; lookup[index3] = (int)(output - outputStart + 3); if (((unsigned char)(next8 >> 0 * 8) == FC_MATCH) && (value > 0)) goto FC_PP_MATCH_FOUND4;

                    prev4 = next4; input += 4; output += 4;

                    continue;
                }

FC_PP_MATCH_FOUND4:
                input += 1; output += 1;
FC_PP_MATCH_FOUND3:
                input += 1; output += 1;
FC_PP_MATCH_FOUND2:
                input += 1; output += 1;
FC_PP_MATCH_FOUND1:
                input += 1;

                if (*input != 255)
                {
                    int len = minLen; while (true) { len += *input; if (*input++ != 254) break; }

                    const unsigned char * reference = outputStart + value;
                          unsigned char * outputEnd = output + len;

                    while (output < outputEnd) { *output++ = *reference++; }

                    prev4 = *(unsigned int *)(output - 4);
                }
                else
                {
                    input++; output++; prev4 = *(unsigned int *)(output - 4); 
                }
            }
        }
#endif

        {
            unsigned int context = output[-1] | (output[-2] << 8) | (output[-3] << 16) | (output[-4] << 24);

            while (input < inputEnd)
            {
                unsigned int index = ((context >> 15) ^ context ^ (context >> 3)) & mask;
                int value = lookup[index]; lookup[index] = (int)(output - outputStart);
                if (*input == FC_MATCH && value > 0)
                {
                    input++;
                    if (*input != 255)
                    {
                        int len = minLen; while (true) { len += *input; if (*input++ != 254) break; }

                        const unsigned char * reference = outputStart + value;
                              unsigned char * outputEnd = output + len;

                        while (output < outputEnd) *output++ = *reference++;

                        context = output[-1] | (output[-2] << 8) | (output[-3] << 16) | (output[-4] << 24);
                    }
                    else
                    {
                        input++; context = (context << 8) | (*output++ = FC_MATCH);
                    }
                }
                else
                {
                    context = (context << 8) | (*output++ = *input++);
                }
            }
        }

        free(lookup);

        return (int)(output - outputStart);
    }

    return FC_LACK_OF_MEMORY;
}



#endif