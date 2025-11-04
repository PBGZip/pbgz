#pragma once

/* 指令集检测*/

#include <stdio.h>
#include <stdint.h>
#include <emmintrin.h>
#include <exception>

#define CPUID_REGION_SHIFT_POS 0
#define CPUID_REGION_SHIFT_REG 10
#define CPUID_REGION_SHIFT_FID_SUB 12
#define CPUID_REGION_SHIFT_FID 20
#define CPUID_REGION_SHIFT_LEN 5

#define CPUID_REGION_MASK_REG 0x00000C00        /* 寄存器: 0:EAX, 1:EBX, 2:ECX, 3:EDX */
#define CPUID_REGION_MASK_FID_SUB 0x000FF000 /* 子功能号(低8位) */
#define CPUID_REGION_MASK_POS 0x0000001F         /* [0,31] bit偏移 */
#define CPUID_REGION_MASK_LEN 0x000003E0         /* [1, 32] bit长 */

#define CPUID_REGION_PRODUCE(f, fs, r, p, l) (((f)&0xF0000000) \
    | ((f) << CPUID_REGION_SHIFT_FID & 0x0FF00000) \
    | ((fs) << CPUID_REGION_SHIFT_FID_SUB & CPUID_REGION_MASK_FID_SUB) \
    | ((r) << CPUID_REGION_SHIFT_REG & CPUID_REGION_MASK_REG) \
    | ((p) << CPUID_REGION_SHIFT_POS & CPUID_REGION_MASK_POS) \
    | (((l)-1) << CPUID_REGION_SHIFT_LEN & CPUID_REGION_MASK_LEN))

#define CPUF_BMI2 CPUID_REGION_PRODUCE(7,0,1,8,1)
#define CPUF_AVX CPUID_REGION_PRODUCE(1, 0, 2, 28, 1)
#define CPUF_AVX2 CPUID_REGION_PRODUCE(7, 0, 1, 5, 1)
#define CPUF_OSXSAVE CPUID_REGION_PRODUCE(1, 0, 2, 27, 1)
#define CPUF_XFSM CPUID_REGION_PRODUCE(0xD, 0, 0, 0, 32)

#define CPUID_REGION_FID(cregion)    ( ((cregion)&0xF0000000) | (((cregion) & 0x0FF00000)>>CPUID_REGION_SHIFT_FID) )
#define CPUID_REGION_FID_SUB(cregion)    ( ((cregion) & CPUID_REGION_MASK_FID_SUB)>>CPUID_REGION_SHIFT_FID_SUB )

class Hardware
{
public:
    Hardware() {}
    static bool isSupportSimd()
    {
        if (getCpuIdRegion(CPUF_BMI2) <= 0)
            return false; /* bmi2 */
        if (getSseLevel() < 2)
            return false; /* sse2 */
        if (getSimdAvxLevel() < 2)
            return false; /* avx2 */
        return true;
    }

private:

    static inline void getCpuId(uint32_t cpuInfo[4], const uint32_t tInfo, const uint32_t vecx)
    {
        asm volatile("cpuid"
                     : "=a"(cpuInfo[0]),
                       "=b"(cpuInfo[1]),
                       "=c"(cpuInfo[2]),
                       "=d"(cpuInfo[3])
                     : "0"(tInfo), "2"(vecx)
                     : "memory");
    }

    /* 根据CPUIDREGION获取CPUID区域 */
    static inline uint32_t getCpuIdRegion(const uint32_t &cInfo)
    {
        uint32_t buffer[4];
        getCpuId(buffer, CPUID_REGION_FID(cInfo), CPUID_REGION_FID_SUB(cInfo));
        return (((buffer[((cInfo & 0x00000C00) >> 10)]) >> ((cInfo & 0x0000001F))) & (((uint32_t)-1) >> (32 - (((cInfo & 0x000003E0) >> 5) + 1))));
    }

    static inline uint32_t getSseLevel()
    {
        uint32_t sseLevel = 0;
        uint32_t cpuInfo[4];

        getCpuId(cpuInfo, 1, 0);
        if (cpuInfo[3] & 0x02000000) /* sse: bit 25 */
        {
            ++sseLevel;
            if (cpuInfo[3] & 0x04000000) /* sse2: bit 26 */
            {
                ++sseLevel;
                if (cpuInfo[2] & 0x00000001) /* sse3: bit 0 */
                {
                    ++sseLevel;
                    if (cpuInfo[2] & 0x00000100) /* ssse3: bit 9 */
                    {
                        ++sseLevel;
                        if (cpuInfo[2] & 0x00080000) /* sse41: bit 19 */
                        {
                            ++sseLevel;
                            if (cpuInfo[2] & 0x00100000) /* sse42: bit 20 */
                                ++sseLevel;
                        }
                    }
                }
            }
        }
        return sseLevel;
    }

    /* 获取AVX系列指令集的支持级别 */
    static inline int getSimdAvxLevel()
    {
        int avxLevel = 0;
        /* 获取处理器信息*/
        if (0 != getCpuIdRegion(CPUF_AVX))
        {
            ++avxLevel;
            if (0 != getCpuIdRegion(CPUF_AVX2))
                ++avxLevel;
        }

        /*  检查操作系统是否支持*/
        if (0 != getCpuIdRegion(CPUF_OSXSAVE)) /* XGETBV enabled for application use */
        {
            uint32_t n = getCpuIdRegion(CPUF_XFSM); /* XCR0: XFeatureSupportedMask register */
            if (6 == (n & 6)) /* XCR0[2:1] =  11b' (XMM state and YMM state are enabled by OS) */
                return avxLevel;
        }
        return 0;
    }
};
