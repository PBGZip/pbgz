/*
 * aux_block_consumer.h - Header file for pbgz project
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

#ifndef _AUX_BLOCK_CONSUMER_H_
#define _AUX_BLOCK_CONSUMER_H_

#include <memory>
#include <vector>

#include "io_block.h"

/*
 * Claimant for auxiliary blocks. The reader asks each registered claimant in turn;
 * each decides whether a block is for it. Unclaimed blocks are skipped. Adding a new
 * auxiliary block type only requires registering another claimant.
 *
 * Threading: claim runs only on the reader thread, but state it writes is read by
 * worker threads, so cross-thread visibility is the claimant's responsibility.
 */
class AuxBlockConsumer {
public:
    virtual ~AuxBlockConsumer() { }

    /*
     * Return true if the block is claimed, false to pass it to the next claimant.
     * packageIndex identifies the pbgz package the block belongs to (the parsed
     * package sequence number, valid for both files and pipes). Claimants must keep
     * blocks per package and not let a later one overwrite an earlier one.
     * Assumes at most one auxiliary block of each kind per package.
     */
    virtual bool claim(RoughIOBlock* blockPtr, int64_t packageIndex) = 0;
};

/* Auxiliary block payload shared by reference: a decoding data block must keep the payload it references alive until it finishes decoding. */
using AuxPayloadPtr = std::shared_ptr<const std::vector<uint8_t> >;

#endif
