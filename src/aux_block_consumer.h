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
 * Claimant for auxiliary blocks.
 *
 * The block reader knows nothing about any concrete auxiliary block type; when it
 * passes an auxiliary block it simply asks each registered claimant in turn, and the
 * claimant itself decides whether the block is meant for it. Adding a new auxiliary
 * block type therefore only requires registering one more claimant - the reader does
 * not change a line. Auxiliary blocks no one claims are silently skipped, which is
 * exactly the forward-compatible behavior old versions need when they encounter the
 * new format.
 *
 * Threading contract: claim is only called on the reader thread, so claimants never
 * run concurrently with each other. However, the state a claimant writes is read by
 * worker threads, so cross-thread visibility is the claimant's own responsibility.
 */
class AuxBlockConsumer {
public:
    virtual ~AuxBlockConsumer() { }

    /*
     * Returning true means this block has been claimed; false means the claimant is
     * not interested and the next claimant should be asked.
     *
     * packageIndex is the sequence number of the pbgz package this block belongs to,
     * i.e. the auxiliary block's identity. Once multiple packages are cat-concatenated,
     * the same auxiliary block type appears more than once, so a claimant must keep
     * them per package and must not let a later one overwrite an earlier one - data
     * blocks are decoded in parallel as they are read, and when a later package's
     * auxiliary block lands, the earlier package's blocks are often still being decoded.
     * Overwriting would silently make them use the wrong data.
     *
     * The "absolute file offset" used to serve as the identity here, and that was
     * wrong: under a non-seekable pipe input the offset collapses to 0, the claimant
     * registers under 0 while data blocks look up by in-package relative offset, and
     * every lookup misses. The package sequence number depends only on the number of
     * parsed package headers, so it holds equally for files and pipes.
     *
     * It is currently assumed that a package contains at most one auxiliary block of
     * each kind, which holds until priors are sharded. Once sharding lands, the key
     * must be extended to (packageIndex, ordinal of this auxiliary block type within
     * the package), and the data block meta must record that ordinal; the offset in the
     * block meta serves only as a seek mechanism and for verification, not as an index.
     */
    virtual bool claim(RoughIOBlock* blockPtr, int64_t packageIndex) = 0;
};

/* Auxiliary block payload shared by reference: a decoding data block must keep the payload it references alive until it finishes decoding. */
using AuxPayloadPtr = std::shared_ptr<const std::vector<uint8_t> >;

#endif
