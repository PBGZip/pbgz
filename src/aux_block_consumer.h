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
 * 辅助块的认领者。
 *
 * 块读者不认识任何具体的辅助块类型，它只负责在路过辅助块时逐个询问已注册的认领者，
 * 由认领者自己判断这块是不是给它的。新增一种辅助块因此只需要多注册一个认领者，
 * 读者一行都不用改；没有人认领的辅助块被安静跳过，这正是老版本读到新格式时
 * 需要的前向兼容行为。
 *
 * 线程约定：claim 只在读线程里被调用，认领者之间不会并发。但认领者写下的状态
 * 会被工作线程读取，所以跨线程可见性由认领者自己负责。
 */
class AuxBlockConsumer {
public:
    virtual ~AuxBlockConsumer() { }

    /*
     * 返回 true 表示本块已被认领；返回 false 表示不感兴趣，继续问下一个认领者。
     *
     * blockAddress 是本块在整个输入流里的绝对起点，也就是它的身份。多个包 cat 到
     * 一起之后同一种辅助块会出现多份，认领者必须按地址分别保存，不能用后来的覆盖
     * 先前的——数据块是边读边并行解的，后一个包的辅助块落地时前一个包的块往往还在
     * 解，覆盖会让它们用错数据且不报错。用绝对地址而不是包起点做键，是因为包起点
     * 只在"每个包恰好一份"时才够用，先验一旦分片就不成立。
     */
    virtual bool claim(RoughIOBlock* blockPtr, int64_t blockAddress) = 0;
};

/* 辅助块载荷按引用共享：解码中的数据块必须能让它引用的那份数据活到自己解完。 */
using AuxPayloadPtr = std::shared_ptr<const std::vector<uint8_t> >;

#endif
