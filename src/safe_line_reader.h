/*
 * safe_line_reader.h - Bounds-safe line iteration over a RoughIOBlock
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

#pragma once

#include <stdint.h>
#include <vector>

#include "io_block.h"

/*
 * SafeLineReader
 *
 * Iterates the newline-delimited lines of a RoughIOBlock in a single place,
 * clamping every offset against the block's data length.  A corrupt or
 * out-of-order npos table can never drive an out-of-bounds read through this
 * API -- malformed lines are silently skipped.
 *
 * This centralises the bounds-safety that field extractors previously
 * duplicated, so the guarantee lives in exactly one location.
 */
class SafeLineReader {
public:
    explicit SafeLineReader(RoughIOBlock* block)
        : buffer_(block ? block->getBuffer() : nullptr),
          npos_(block ? &block->getNpos() : nullptr),
          dataLen_(block ? (uint64_t)block->getDataLen() : 0),
          lineIdx_(0),
          scanned_(0)
    {
        if (npos_ == nullptr || buffer_ == nullptr || dataLen_ == 0) {
            lineCount_ = 0;
        } else {
            lineCount_ = (uint32_t)npos_->size();
        }
    }

    /*
     * Advance to the next line.  On success returns true and sets line/len
     * to a view of the line content (excluding the trailing newline).
     * Returns false when lines are exhausted.
     */
    bool nextLine(const uint8_t*& line, uint32_t& len)
    {
        while (lineIdx_ < lineCount_) {
            uint64_t start = (lineIdx_ == 0) ? 0 : (uint64_t)(*npos_)[lineIdx_ - 1] + 1;
            uint64_t end = (*npos_)[lineIdx_];
            ++lineIdx_;

            if (end > dataLen_) {
                end = dataLen_;
            }
            if (start >= end) {
                continue;
            }
            line = buffer_ + start;
            len = (uint32_t)(end - start);
            scanned_ += end - start + 1;
            return true;
        }
        return false;
    }

    uint64_t scannedBytes() const { return scanned_; }

private:
    const uint8_t* buffer_;
    const std::vector<uint32_t>* npos_;
    uint64_t dataLen_;
    uint32_t lineIdx_;
    uint32_t lineCount_;
    uint64_t scanned_;
};
