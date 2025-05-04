#ifndef _ACTG_H_
#define _ACTG_H_

#include <stdint.h>
#include <mutex>
#include <bitset>
#include "manager.h"

/* ACTG pair对应的表 */
static uint16_t *actgp = nullptr;
static std::mutex actgp_mutex;

/*  ACTG -> 00011011*/
static inline int64_t actg_squash(const uint8_t *base, int64_t base_len, uint8_t *dst)
{
    uint8_t *squash;
    int64_t n, len4align, len4remain;

    squash = dst;
    len4align = (base_len >> 2) << 2;
    len4remain = base_len - len4align;

    for (n = 0; n < len4align; n += 4)
        *squash++ = ((base[n + 3] & 0x06) >> 1) | ((base[n + 2] & 0x06) << 1) | ((base[n + 1] & 0x06) << 3) | ((base[n] & 0x06) << 5);

    if (len4remain > 0)
    {
        *squash = 0;
        for (n = 0; n < len4remain; n++)
            *squash |= ((base[len4align + n] & 0x06) >> 1) << (6 - (n << 1));
        squash++;
    }

    return (squash - dst);
}

/* 求碱基的互补链 */
static inline void actg_pair(uint8_t *dst, const uint8_t *src, const size_t len)
{
    if (!actgp) /* 为空时创建加速查找表 */
    {
        std::unique_lock<std::mutex> guard(actgp_mutex);
        if (!actgp)
        {
            safe_alloc(65536, uint16_t, actgp);
            uint32_t i1, i2, i3, i4, k;
            uint8_t v[4] = {0, 2, 4, 6};
            uint8_t actg[5];
            static uint8_t actgr7[7] = {'T', 'T', 'G', 'G', 'A', 'A', 'C'};
            actg[4] = '\0';

            for (i1 = 0; i1 < 4; i1++)
            {
                for (i2 = 0; i2 < 4; i2++)
                {
                    k = (v[i1] << 8) | v[i2];
                    sprintf((char *)actg, "%c%c", actgr7[v[i1]], actgr7[v[i2]]);
                    actgp[k] = *((uint16_t *)(actg));
                }
            }
        }
    }
    uint8_t *ps, *pd;
    uint64_t x;
    uint8_t *ealign8 = (uint8_t *)src + (len >> 3 << 3);
    uint8_t *e = (uint8_t *)src + len - 1;
    pd = dst + len;
    for (ps = (uint8_t *)src; ps != ealign8; ps += 8)
    {
        x = (*((uint64_t *)ps)) & 0x0606060606060606;

        pd -= 2;
        *((uint16_t *)(pd)) = actgp[*((uint16_t *)(&x))];
        pd -= 2;
        *((uint16_t *)(pd)) = actgp[*((uint16_t *)(&x) + 1)];
        pd -= 2;
        *((uint16_t *)(pd)) = actgp[*((uint16_t *)(&x) + 2)];
        pd -= 2;
        *((uint16_t *)(pd)) = actgp[*((uint16_t *)(&x) + 3)];
    }

    ps = ealign8;
    if (ps <= e)
    {
        uint8_t actgr7[7] = {'T', 'T', 'G', 'G', 'A', 'A', 'C'};
        do
        {
            *--pd = actgr7[*ps & 0x6];
            ps++;
        } while (ps <= e);
    }
}

/* 求两个字符串不同字符的个数 */
static inline uint32_t get_diff_cnt(const uint8_t *s1, const uint8_t *s2, uint32_t len)
{
    uint32_t n = 0, zcnt = 0;
    for (n = 0; n < len; n++)
        zcnt += !((*(s1 + n)) ^ (*(s2 + n)));
    return len - zcnt;
}

/* 求squash后的actg不同碱基的个数，如01001011与01101110的个数为3 */
static inline uint32_t actg_squash_diffcnt(const uint8_t *s1, const uint8_t *s2, uint32_t len)
{
    uint8_t x8;
    uint64_t x64;
    uint32_t cnt = 0, n, next;
    uint32_t align8 = len >> 3 << 3;
    for (n = 0, next = 8; next < align8; n += 8)
    {
        x64 = (*((uint64_t *)(s1 + n))) ^ (*((uint64_t *)(s2 + n)));
        x64 = (x64 & 0x5555555555555555) | ((x64 >> 1) & 0x5555555555555555);
        cnt += std::bitset<64>{x64}.count();
        next += 8;
    }

    if (n == len)
        return cnt;

    for (; n < len; n++)
    {
        x8 = (*(s1 + n)) ^ (*(s2 + n));
        x8 = (x8 & 0x55) | ((x8 >> 1) & 0x55);
        cnt += std::bitset<8>{x8}.count();
    }
    return cnt;
}

/* 计算squash之后的base与reference的mapping关系，2bits放到一个字节的末尾，如果内容相同则存0，否则存base原始的2bits */
static inline uint32_t actg_stretch_mapping(const uint8_t *squash_base, const uint8_t *squash_refe, uint32_t squash_len, uint8_t *dst)
{
    uint32_t n, m, offset = 0, next;
    uint8_t chb, chr;

    for (n = 0, next = 2; next <= squash_len; n += 2)
    {
        if ((*((uint16_t *)(squash_base + n))) != *((uint16_t *)(squash_refe + n)))
        {
            /* first byte */
            if ((*(squash_base + n)) == (*(squash_refe + n)))
            {
                *((uint32_t *)(dst + offset)) = 0;
                offset += 4;
            }
            else
            {
                chb = *(squash_base + n);
                chr = *(squash_refe + n);
                *(dst + offset++) = ((chb & 0xC0) == (chr & 0xC0)) ? 0 : ((chb & 0xC0) >> 6);
                *(dst + offset++) = ((chb & 0x30) == (chr & 0x30)) ? 0 : ((chb & 0x30) >> 4);
                *(dst + offset++) = ((chb & 0xC) == (chr & 0xC)) ? 0 : ((chb & 0xC) >> 2);
                *(dst + offset++) = ((chb & 0x3) == (chr & 0x3)) ? 0 : (chb & 0x3);
            }
            /* second byte */
            if ((*(squash_base + n + 1)) == (*(squash_refe + n + 1)))
            {
                *((uint32_t *)(dst + offset)) = 0;
                offset += 4;
            }
            else
            {
                chb = *(squash_base + n + 1);
                chr = *(squash_refe + n + 1);
                *(dst + offset++) = ((chb & 0xC0) == (chr & 0xC0)) ? 0 : ((chb & 0xC0) >> 6);
                *(dst + offset++) = ((chb & 0x30) == (chr & 0x30)) ? 0 : ((chb & 0x30) >> 4);
                *(dst + offset++) = ((chb & 0xC) == (chr & 0xC)) ? 0 : ((chb & 0xC) >> 2);
                *(dst + offset++) = ((chb & 0x3) == (chr & 0x3)) ? 0 : (chb & 0x3);
            }
        }
        else
        {
            *((uint64_t *)(dst + offset)) = 0;
            offset += 8;
        }
        next += 2;
    }
    if (n == squash_len)
        return offset;

    /* last byte*/
    if ((*(squash_base + n)) == (*(squash_refe + n)))
    {
        *((uint32_t *)(dst + offset)) = 0;
        offset += 4;
    }
    else
    {
        chb = *(squash_base + n);
        chr = *(squash_refe + n);
        *(dst + offset++) = ((chb & 0xC0) == (chr & 0xC0)) ? 0 : ((chb & 0xC0) >> 6);
        *(dst + offset++) = ((chb & 0x30) == (chr & 0x30)) ? 0 : ((chb & 0x30) >> 4);
        *(dst + offset++) = ((chb & 0xC) == (chr & 0xC)) ? 0 : ((chb & 0xC) >> 2);
        *(dst + offset++) = ((chb & 0x3) == (chr & 0x3)) ? 0 : (chb & 0x3);
    }
    return offset;
}

/* 计算squash之后的base与reference的mapping关系，2bits放到一个字节的末尾 */
static inline uint32_t actg_stretch_mapping_xor(const uint8_t *squash_base, const uint8_t *squash_refe, uint32_t squash_len, uint8_t *dst)
{
    uint32_t n, m, offset = 0, x32, next;
    uint64_t x64, x64_1, x64_2, p64;
    uint8_t *p1, *p2, x, *p;
    uint16_t x16;

#define STRETCH_XOR_2BYTE(n)                                                                                                             \
    {                                                                                                                                    \
        x16 = (*((uint16_t *)(squash_base + n))) ^ (*((uint16_t *)(squash_refe + n)));                                                   \
        p = (uint8_t *)(&x16);                                                                                                           \
        x64_1 = (uint64_t)(*p++);                                                                                                        \
        x64_2 = (uint64_t)(*p++);                                                                                                        \
        *((uint64_t *)(dst + offset)) = ((x64_1 & 0xC0) >> 6) | ((x64_1 & 0x30) << 4) | ((x64_1 & 0xC) << 14) | ((x64_1 & 0x3) << 24) |  \
                                        ((x64_2 & 0xC0) << 26) | ((x64_2 & 0x30) << 36) | ((x64_2 & 0xC) << 46) | ((x64_2 & 0x3) << 56); \
        offset += 8;                                                                                                                     \
    }

#define STRETCH_XOR_1BYTE(n)                                                                                                   \
    {                                                                                                                          \
        x32 = (uint32_t)((*(squash_base + n)) ^ (*(squash_refe + n)));                                                         \
        *((uint32_t *)(dst + offset)) = ((x32 & 0xC0) >> 6) | ((x32 & 0x30) << 4) | ((x32 & 0xC) << 14) | ((x32 & 0x3) << 24); \
        offset += 4;                                                                                                           \
    }

    for (n = 0, next = 8; next <= squash_len; n += 8)
    {
        x64 = (*((uint64_t *)(squash_base + n))) ^ (*((uint64_t *)(squash_refe + n)));
        p = (uint8_t *)(&x64);

        x64_1 = (uint64_t)(*p++);
        x64_2 = (uint64_t)(*p++);
        *((uint64_t *)(dst + offset)) = ((x64_1 & 0xC0) >> 6) | ((x64_1 & 0x30) << 4) | ((x64_1 & 0xC) << 14) | ((x64_1 & 0x3) << 24) |
                                        ((x64_2 & 0xC0) << 26) | ((x64_2 & 0x30) << 36) | ((x64_2 & 0xC) << 46) | ((x64_2 & 0x3) << 56);
        offset += 8;

        x64_1 = (uint64_t)(*p++);
        x64_2 = (uint64_t)(*p++);
        *((uint64_t *)(dst + offset)) = ((x64_1 & 0xC0) >> 6) | ((x64_1 & 0x30) << 4) | ((x64_1 & 0xC) << 14) | ((x64_1 & 0x3) << 24) |
                                        ((x64_2 & 0xC0) << 26) | ((x64_2 & 0x30) << 36) | ((x64_2 & 0xC) << 46) | ((x64_2 & 0x3) << 56);
        offset += 8;

        x64_1 = (uint64_t)(*p++);
        x64_2 = (uint64_t)(*p++);
        *((uint64_t *)(dst + offset)) = ((x64_1 & 0xC0) >> 6) | ((x64_1 & 0x30) << 4) | ((x64_1 & 0xC) << 14) | ((x64_1 & 0x3) << 24) |
                                        ((x64_2 & 0xC0) << 26) | ((x64_2 & 0x30) << 36) | ((x64_2 & 0xC) << 46) | ((x64_2 & 0x3) << 56);
        offset += 8;

        x64_1 = (uint64_t)(*p++);
        x64_2 = (uint64_t)(*p++);
        *((uint64_t *)(dst + offset)) = ((x64_1 & 0xC0) >> 6) | ((x64_1 & 0x30) << 4) | ((x64_1 & 0xC) << 14) | ((x64_1 & 0x3) << 24) |
                                        ((x64_2 & 0xC0) << 26) | ((x64_2 & 0x30) << 36) | ((x64_2 & 0xC) << 46) | ((x64_2 & 0x3) << 56);
        offset += 8;

        next += 8;
    }

    switch (squash_len - n)
    {
    case 0:
        break;
    case 1:
        STRETCH_XOR_1BYTE(n);
        break;
    case 2:
        STRETCH_XOR_2BYTE(n);
        break;
    case 3:
        STRETCH_XOR_2BYTE(n);
        n += 2;
        STRETCH_XOR_1BYTE(n);
        break;
    case 4:
        STRETCH_XOR_2BYTE(n);
        n += 2;
        STRETCH_XOR_2BYTE(n);
        break;
    case 5:
        STRETCH_XOR_2BYTE(n);
        n += 2;
        STRETCH_XOR_2BYTE(n);
        n += 2;
        STRETCH_XOR_1BYTE(n);
        break;
    case 6:
        STRETCH_XOR_2BYTE(n);
        n += 2;
        STRETCH_XOR_2BYTE(n);
        n += 2;
        STRETCH_XOR_2BYTE(n);
        break;
    case 7:
        STRETCH_XOR_2BYTE(n);
        n += 2;
        STRETCH_XOR_2BYTE(n);
        n += 2;
        STRETCH_XOR_2BYTE(n);
        n += 2;
        STRETCH_XOR_1BYTE(n);
        break;
    default:
        break;
    }
    return offset;
}

/* 将ACTG编码，编码后一个碱基还是占1个字节，但是只有最后两个bits有效 */
static inline void actg_encode(const uint8_t *src, uint8_t *dst, uint32_t len)
{
    uint32_t n, align8 = len >> 3 << 3;

    for (n = 0; n < align8; n += 8)
        *((uint64_t *)(dst + n)) = ((*((uint64_t *)(src + n))) >> 1) & 0x303030303030303;

    if (n == len)
        return;

    for (; n < len; n++)
        *(dst + n) = ((*(src + n)) >> 1) & 0x3;
}

/* 以字节为单位异或 */
static inline void actg_xor(const uint8_t *x1, const uint8_t *x2, uint8_t *out, uint32_t len)
{
    uint32_t n, align8 = (len >> 3) << 3;
    for (n = 0; n < align8; n += 8)
        *((uint64_t *)(out + n)) = (*((uint64_t *)(x1 + n))) ^ (*((uint64_t *)(x2 + n)));

    switch (len - n)
    {
    case 0:
        break;
    case 1:
        *(out + n) = (*(x1 + n)) ^ (*(x2 + n));
        break;
    case 2:
        *((uint16_t *)(out + n)) = (*((uint16_t *)(x1 + n))) ^ (*((uint16_t *)(x2 + n)));
        break;
    case 3:
        *((uint16_t *)(out + n)) = (*((uint16_t *)(x1 + n))) ^ (*((uint16_t *)(x2 + n)));
        n += 2;
        *(out + n) = (*(x1 + n)) ^ (*(x2 + n));
        break;
    case 4:
        *((uint32_t *)(out + n)) = (*((uint32_t *)(x1 + n))) ^ (*((uint32_t *)(x2 + n)));
        break;
    case 5:
        *((uint32_t *)(out + n)) = (*((uint32_t *)(x1 + n))) ^ (*((uint32_t *)(x2 + n)));
        n += 4;
        *(out + n) = (*(x1 + n)) ^ (*(x2 + n));
        break;
    case 6:
        *((uint32_t *)(out + n)) = (*((uint32_t *)(x1 + n))) ^ (*((uint32_t *)(x2 + n)));
        n += 4;
        *((uint16_t *)(out + n)) = (*((uint16_t *)(x1 + n))) ^ (*((uint16_t *)(x2 + n)));
        break;
    case 7:
        *((uint32_t *)(out + n)) = (*((uint32_t *)(x1 + n))) ^ (*((uint32_t *)(x2 + n)));
        n += 4;
        *((uint16_t *)(out + n)) = (*((uint16_t *)(x1 + n))) ^ (*((uint16_t *)(x2 + n)));
        n += 2;
        *(out + n) = (*(x1 + n)) ^ (*(x2 + n));
        break;
    default:
        break;
    }
}

#endif
