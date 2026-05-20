/*
 * codec_actuator_adapter.h - Head file for pbgz project
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

class PbgzEngine;

#include "actuator.h"
#include "binary_actuator.h"
#include "fastq_actuator.h"
#include "sam_actuator.h"
#include "utils/memory_util.h"

class BinaryCompressActuator : public Actuator {
public:
    BinaryCompressActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr) : Actuator(inPtr, outPtr, engine) {

    }

    virtual ~BinaryCompressActuator() {
        MemoryUtil::safeDeleteClass(codecActuator);
    }

    int32_t process() {
        if (codecActuator != nullptr) {
            return codecActuator->compress();
        }

        return -1;
    }

    int32_t initial() {
        codecActuator = MemoryUtil::safeNewClass<BinaryCodecActuator>(inBlockPtr, outBlockPtr, pbgzEngine);
        if (codecActuator == nullptr) {
            return -1;
        }
        return 0;
    }

    BinaryCodecActuator* getCodecActuator() {
        return codecActuator;
    }

protected:
    BinaryCodecActuator* codecActuator;
};


class BinaryDecompressActuator : public Actuator {
public:
    BinaryDecompressActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine = nullptr) : Actuator(inPtr, outPtr, engine) {

    }

    virtual ~BinaryDecompressActuator() {
        MemoryUtil::safeDeleteClass(codecActuator);
    }

    int32_t process() {
        if (codecActuator != nullptr) {
            return codecActuator->decompress();
        }

        return -1;
    }

    int32_t initial() {
        codecActuator = MemoryUtil::safeNewClass<BinaryCodecActuator>(inBlockPtr, outBlockPtr, pbgzEngine);
        if (codecActuator == nullptr) {
            return -1;
        }
        return 0;
    }

    BinaryCodecActuator* getCodecActuator() {
        return codecActuator;
    }

protected:
    BinaryCodecActuator* codecActuator;
};

template <typename T>
class CompressActuator : public Actuator {
public:
    CompressActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine, Reference* pRef) : Actuator(inPtr, outPtr, engine) {
        pReference = pRef;
    }

    virtual ~CompressActuator() {
        MemoryUtil::safeDeleteClass(codecActuator);
    }

    int32_t process() {
        if (codecActuator != nullptr) {
            return codecActuator->compress();
        }

        return -1;
    }

    int32_t initial() {
        codecActuator = MemoryUtil::safeNewClass<T>(inBlockPtr, outBlockPtr, pbgzEngine, pReference);
        if (codecActuator == nullptr) {
            return -1;
        }

        return 0;
    }

    T* getCodecActuator() {
        return codecActuator;
    }

protected:
    T* codecActuator;
    Reference* pReference;
};


template <typename T>
class DecompressActuator : public Actuator {
public:
    DecompressActuator(RoughIOBlock* inPtr, RoughIOBlock* outPtr, PbgzEngine* engine, Reference* pRef) : Actuator(inPtr, outPtr, engine) {
        pReference = pRef;
    }

    virtual ~DecompressActuator() {
        MemoryUtil::safeDeleteClass(codecActuator);
    }

    int32_t process() {
        if (codecActuator != nullptr) {
            return codecActuator->decompress();
        }

        return -1;
    }

    int32_t initial() {
        codecActuator = MemoryUtil::safeNewClass<T>(inBlockPtr, outBlockPtr, pbgzEngine, pReference);
        if (codecActuator == nullptr) {
            return -1;
        }
        return 0;
    }

    T* getCodecActuator() {
        return codecActuator;
    }

protected:
    T* codecActuator;
    Reference* pReference;
};

using FastqCompressActuator = CompressActuator<FastqCodecActuator>;
using FastqDecompressActuator = DecompressActuator<FastqCodecActuator>;


using SamCompressActuator = CompressActuator<SamCodecActuator>;
using SamDecompressActuator = DecompressActuator<SamCodecActuator>;

