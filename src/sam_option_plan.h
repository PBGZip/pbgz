/*
 * sam_option_plan.h - file-level SAM OPTION layout plan
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

#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

class Reference;

/* Build the per-tag layout plan from one block (runs on a reader-side thread).  `mode` is the
   compression profile (-m), so the trial encodes exactly as the blocks will.  For every tag it
   records the layout that won the exhaustive trial.  Returns the number of tags planned; 0
   means no plan (blocks then trial every layout). */
int32_t buildSamOptionPlan(uint8_t* buffer, uint32_t dataLen,
                           const std::vector<size_t>& npos, Reference* ref, uint8_t mode,
                           std::unordered_map<std::string, std::string>& out);

/*
 * Per-file OPTION layout plan, built once from the first data block.  It records the layout
 * that won the exhaustive trial for every tag; blocks then encode that layout directly and
 * skip the per-block comparison.  A tag the plan does not cover (or a block where the planned
 * layout cannot be produced) falls back to the full trial.
 *
 * The plan is built before the coder workers need it and is read-only afterwards, so the
 * archive does not depend on thread scheduling.
 */
struct SamOptionPlan {
    std::mutex mutex;
    std::condition_variable cv;
    /* Set (under the mutex) before the block that triggered the build is dispatched, so every
       worker that sees the plan knows it will become ready. */
    bool started = false;
    bool ready = false;
    /* Key: tag name + '\x1f' + tag type.  Value: the layout to write. */
    std::unordered_map<std::string, std::string> layout;
};
