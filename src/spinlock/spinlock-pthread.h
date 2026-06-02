/*
 * spinlock-pthread.h - Spinlock using pthread mutex
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

#ifndef _SPINLOCK_PTHREAD_H
#define _SPINLOCK_PTHREAD_H

#define SPINLOCK_ATTR static __inline __attribute__((always_inline, no_instrument_function))

#define spinlock pthread_mutex_t

SPINLOCK_ATTR void spin_lock(spinlock *lock)
{
    pthread_mutex_lock(lock);
}

SPINLOCK_ATTR void spin_unlock(spinlock *lock)
{
    pthread_mutex_unlock(lock);
}

#define SPINLOCK_INITIALIZER { 0 }

#endif /* _SPINLOCK_H */
