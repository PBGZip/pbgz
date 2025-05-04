#ifndef _HARDWARE_H_
#define _HARDWARE_H_

/* 指令集检测*/

#include <stdio.h>
#include <stdint.h>
#include <emmintrin.h>
#include <exception>

#define CPUIDREGION_SHIFT_POS 0
#define CPUIDREGION_SHIFT_REG 10
#define CPUIDREGION_SHIFT_FIDSUB 12
#define CPUIDREGION_SHIFT_FID 20
#define CPUIDREGION_SHIFT_LEN 5

#define CPUIDREGION_MASK_REG 0x00000C00        /* 寄存器: 0:EAX, 1:EBX, 2:ECX, 3:EDX */
#define CPUIDREGION_MASK_FIDSUB 0x000FF000 /* 子功能号(低8位) */
#define CPUIDREGION_MASK_POS 0x0000001F         /* [0,31] bit偏移 */
#define CPUIDREGION_MASK_LEN 0x000003E0         /* [1, 32] bit长 */

#define CPUIDREGION_PRODUCE(f, fs, r, p, l) (((f)&0xF0000000) \
    | ((f) << CPUIDREGION_SHIFT_FID & 0x0FF00000) \
    | ((fs) << CPUIDREGION_SHIFT_FIDSUB & CPUIDREGION_MASK_FIDSUB) \
    | ((r) << CPUIDREGION_SHIFT_REG & CPUIDREGION_MASK_REG) \
    | ((p) << CPUIDREGION_SHIFT_POS & CPUIDREGION_MASK_POS) \
    | (((l)-1) << CPUIDREGION_SHIFT_LEN & CPUIDREGION_MASK_LEN))

#define CPUF_BMI2 CPUIDREGION_PRODUCE(7,0,1,8,1)
#define CPUF_AVX CPUIDREGION_PRODUCE(1, 0, 2, 28, 1)
#define CPUF_AVX2 CPUIDREGION_PRODUCE(7, 0, 1, 5, 1)
#define CPUF_OSXSAVE CPUIDREGION_PRODUCE(1, 0, 2, 27, 1)
#define CPUF_XFSM CPUIDREGION_PRODUCE(0xD, 0, 0, 0, 32)

#define CPUIDREGION_FID(cregion)    ( ((cregion)&0xF0000000) | (((cregion) & 0x0FF00000)>>CPUIDREGION_SHIFT_FID) )
#define CPUIDREGION_FIDSUB(cregion)    ( ((cregion) & CPUIDREGION_MASK_FIDSUB)>>CPUIDREGION_SHIFT_FIDSUB )

class hardware
{
public:
    hardware() {}
    bool support_simd()
    {
        if (get_cpuidregion(CPUF_BMI2) <= 0)
            return false; /* bmi2 */
        if (get_sselevel() < 2)
            return false; /* sse2 */
        // if (get_simdavx_level(&nhwavx) < 1)
        //     return false; /* avx */
        if (get_simdavx_level() < 2)
            return false; /* avx2 */
        return true;
    }

private:

    inline void __get_cpuid(uint32_t CInfo[4], const uint32_t TInfo, const uint32_t VECX)
    {
        asm volatile("cpuid"
                     : "=a"(CInfo[0]),
                       "=b"(CInfo[1]),
                       "=c"(CInfo[2]),
                       "=d"(CInfo[3])
                     : "0"(TInfo), "2"(VECX)
                     : "memory");
    }

    /* 根据CPUIDREGION获取CPUID区域 */
    inline uint32_t get_cpuidregion(const uint32_t &CInfo)
    {
        uint32_t buffer[4];
        __get_cpuid(buffer, CPUIDREGION_FID(CInfo), CPUIDREGION_FIDSUB(CInfo));
        return (((buffer[((CInfo & 0x00000C00) >> 10)]) >> ((CInfo & 0x0000001F))) & (((uint32_t)-1) >> (32 - (((CInfo & 0x000003E0) >> 5) + 1))));
    }

    inline uint32_t get_sselevel()
    {
        uint32_t sse_level = 0;
        uint32_t CInfo[4];

        __get_cpuid(CInfo, 1, 0);
        if (CInfo[3] & 0x02000000) /* sse: bit 25 */
        {
            ++sse_level;
            if (CInfo[3] & 0x04000000) /* sse2: bit 26 */
            {
                ++sse_level;
                if (CInfo[2] & 0x00000001) /* sse3: bit 0 */
                {
                    ++sse_level;
                    if (CInfo[2] & 0x00000100) /* ssse3: bit 9 */
                    {
                        ++sse_level;
                        if (CInfo[2] & 0x00080000) /* sse41: bit 19 */
                        {
                            ++sse_level;
                            if (CInfo[2] & 0x00100000) /* sse42: bit 20 */
                                ++sse_level;
                        }
                    }
                }
            }
        }
        return sse_level;
    }

    /* 获取AVX系列指令集的支持级别 */
    inline int get_simdavx_level()
    {
        int avx_level = 0;
        /* 获取处理器信息*/
        if (0 != get_cpuidregion(CPUF_AVX))
        {
            ++avx_level;
            if (0 != get_cpuidregion(CPUF_AVX2))
                ++avx_level;
        }

        /*  检查操作系统是否支持*/
        if (0 != get_cpuidregion(CPUF_OSXSAVE)) /* XGETBV enabled for application use */
        {
            uint32_t n = get_cpuidregion(CPUF_XFSM); /* XCR0: XFeatureSupportedMask register */
            if (6 == (n & 6)) /* XCR0[2:1] =  11b' (XMM state and YMM state are enabled by OS) */
                return avx_level;
        }
        return 0;
    }
};

#endif