# PBGZ 并行化架构设计

## 整体架构图

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                                    PBGZ 并行压缩系统                              │
└─────────────────────────────────────────────────────────────────────────────────┘

┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   输入文件      │    │   线程池管理     │    │   输出文件      │
│  (File/Stdin)  │    │  (CPU核心数)     │    │ (File/Stdout)   │
└─────────────────┘    └─────────────────┘    └─────────────────┘
         │                       │                       │
         ▼                       │                       ▼
┌─────────────────┐              │              ┌─────────────────┐
│   读取线程      │              │              │   写入线程      │
│  (ReadTask)     │◄─────────────┼─────────────►│  (WriteTask)    │
│   单线程        │              │              │    单线程       │
└─────────────────┘              │              └─────────────────┘
         │                       │                       │
         ▼                       │                       ▼
┌─────────────────┐              │              ┌─────────────────┐
│  freeInputPool  │              │              │ outputDataPool  │
│   (空闲输入块)  │              │              │   (输出数据)    │
│   容量: N       │              │              │   容量: 2N      │
└─────────────────┘              │              └─────────────────┘
         │▲                      │                       │▲
          │                      │                        │
         ▼│                      │                        ││
┌─────────────────┐              │              ┌─────────────────┐
│  inputDataPool  │              │              │ freeOutputPool  │
│   (待处理数据)  │              │              │   (空闲输出块)  │
│   容量: N       │              │              │   容量: 2N      │
└─────────────────┘              │              └─────────────────┘
         │▲                      │                       │▲
          │                      │                        │
         ▼│                      │                        ││
    ┌─────────┼─────────┐        │        ┌─────────┼─────────┐
    │         │         │        │        │         │         │
    ▼         ▼         ▼        │        ▼         ▼         ▼
┌─────────┐ ┌─────────┐ ┌─────────┐ │ ┌─────────┐ ┌─────────┐ ┌─────────┐
│Coder    │ │Coder    │ │Coder    │ │ │Coder    │ │Coder    │ │Coder    │
│Thread-0 │ │Thread-1 │ │Thread-N │ │ │Thread-0 │ │Thread-1 │ │Thread-N │
│(压缩线程)│ │(压缩线程)│ │(压缩线程)│ │ │(压缩线程)│ │(压缩线程)│ │(压缩线程)│
└─────────┘ └─────────┘ └─────────┘ │ └─────────┘ └─────────┘ └─────────┘
    │         │         │        │    │         │         │         │
    └─────────┼─────────┘        │    └─────────┼─────────┘         │
              │                  │              │                   │
              ▼                  │              ▼                   ▼
        ┌─────────────┐          │        ┌─────────────┐     ┌─────────────┐
        │   压缩器     │          │        │   排序缓存   │     │   块写入器   │
        │  (FC/BWT)   │          │        │ (顺序保证)   │     │ (PBGZ格式)   │
        └─────────────┘          │        └─────────────┘     └─────────────┘
                                   │
                           ┌─────────────┐
                           │  队列管理器   │
                           │ (阻塞队列)   │
                           └─────────────┘
```

## 核心组件设计

### 1. 线程模型
```
主线程 (main) - 执行startReadTask()，被命名为"readtask"
├── 压缩线程池 (CoderTask) - N个 (N=CPU核心数)
└── 写入线程 (WriteTask) - 1个

总计: N+2个线程
- 1个主线程 (在startReadTask()中重命名为"readtask")
- N个压缩线程 (在startCoderTask()中创建，命名为"codertask_0"到"codertask_N-1")  
- 1个写入线程 (在startWriteTask()中创建，命名为"writetask")

执行顺序:
1. startWriteTask() - 创建写入线程
2. startCoderTask() - 创建N个压缩线程  
3. startReadTask() - 在主线程执行读取任务
4. 等待压缩线程完成
5. 通知写入线程结束
6. 等待写入线程完成
```

### 2. 队列系统
```
┌─────────────────────────────────────────────────────────────┐
│                    阻塞队列网络                              │
├─────────────────┬─────────────────┬─────────────────────────┤
│  freeInputPool  │  inputDataPool  │  freeOutputPool         │
│   (空闲输入块)  │   (待处理数据)  │   (空闲输出块)          │
│   容量: N       │   容量: N       │   容量: 2N              │
└─────────────────┴─────────────────┴─────────────────────────┤
│                  outputDataPool                            │
│                (输出数据队列)                               │
│                  容量: 2N                                  │
└─────────────────────────────────────────────────────────────┘
```

### 3. 数据流向
```
输入数据 → [读取线程] → inputDataPool → [压缩线程池] → outputDataPool → [写入线程] → 输出数据
    ↑                              ↓                              ↑
freeInputPool ← [压缩完成] ← freeOutputPool ← [写入完成] ← [排序缓存]
```

## 并行化设计逻辑

### 1. 生产者-消费者模型

#### **三级流水线设计**:
1. **读取阶段** (Producer): 单线程读取数据块
2. **处理阶段** (Consumer): 多线程并行压缩
3. **写入阶段** (Consumer): 单线程顺序写入

#### **队列容量设计**:
```cpp
// 队列容量与线程数相关
freeInputPool.setCapility(parameter.threadNum);      // N
inputDataPool.setCapility(parameter.threadNum);      // N  
freeOutputPool.setCapility(parameter.threadNum << 1); // 2N
outputDataPool.setCapility(parameter.threadNum << 1); // 2N
```

### 2. 负载均衡策略

#### **工作窃取模式**:
```cpp
// 每个压缩线程独立从队列获取任务
RoughIOBlock* inBlockPtr = inputDataPool.get();  // 阻塞获取
```

#### **动态任务分配**:
- 压缩线程竞争获取数据块
- 快的线程自动处理更多数据
- 避免线程间负载不均

### 3. 同步机制

#### **阻塞队列**:
```cpp
template<typename T>
class BlockingQueue {
    void push(const T &item) {
        std::unique_lock<std::mutex> lock(mutex);
        conditonVar.wait(lock, [this]{return dataQueue.size() < maxQueueSize;});
        dataQueue.push(std::move(item));
        conditonVar.notify_one();
    }
    
    T get() {
        std::unique_lock<std::mutex> lock(mutex);
        conditonVar.wait(lock, [this]{return !this->dataQueue.empty();});
        T value = std::move(dataQueue.front());
        dataQueue.pop();
        return value;
    }
};
```

#### **线程同步**:
- **条件变量**: 队列空/满时自动阻塞
- **互斥锁**: 保护队列操作
- **原子操作**: 计数器和状态管理

### 4. 内存管理

#### **对象池模式**:
```cpp
// 预分配数据块，避免频繁内存分配
for (uint32_t i = 0; i < freeInputPool.getCapility(); ++i) {
    RoughIOBlock* inPtr = new RoughIOBlock();  // 256MB块
    freeInputPool.push(inPtr);
}
```

#### **循环利用**:
- 数据块在队列间循环流动
- 压缩完成后归还到空闲队列
- 减少内存分配/释放开销

### 5. 顺序保证机制

#### **块ID管理**:
```cpp
// 读取时分配递增ID
blockPtr->setBlockId(blockId++);

// 写入时按ID排序
outputSortedCache.sort([](const RoughIOBlock* p1, RoughIOBlock* p2) {
    return p1->getBlockId() <= p2->getBlockId();
});
```

#### **乱序处理**:
- 压缩线程并行处理，输出乱序
- 写入线程缓存并排序
- 保证最终输出顺序正确

### 6. 错误处理和优雅关闭

#### **结束信号传播**:
```cpp
// 读取结束: 向所有压缩线程发送空指针
for (int i = 0; i < parameter.threadNum; ++i) {
    inputDataPool.push(nullptr);  // 结束信号
}

// 压缩结束: 向写入线程发送空指针  
outputDataPool.push(nullptr);    // 结束信号
```

#### **资源清理**:
- 线程join等待完成
- 队列元素逐个清理
- 内存对象正确释放

## 性能优化特点

### 1. **零拷贝设计**
- 数据块在内存中传递
- 避免不必要的数据复制

### 2. **缓存友好**
- 大块数据(256MB)提高缓存命中率
- 减少内存碎片

### 3. **CPU利用率最大化**
- 线程数=CPU核心数
- 避免线程切换开销

### 4. **内存使用可控**
- 队列容量限制内存使用
- 对象池避免内存泄漏

## 扩展性设计

### 1. **可配置线程数**
```cpp
threadNum = sysconf(_SC_NPROCESSORS_ONLN);  // 自动检测
// 或通过 -t 参数手动指定
```

### 2. **可扩展压缩器**
- 插件式压缩器架构
- 支持不同数据类型优化

### 3. **队列容量可调**
- 根据系统内存调整队列大小
- 平衡内存使用和性能

这个并行化设计充分利用了现代多核CPU的计算能力，通过精心设计的队列系统和线程模型，实现了高效的并行压缩处理。
