# PBGZ VTune Profiling Tool

这个工具包提供了使用 Intel VTune Profiler 分析 pbgz 性能的完整解决方案。

## 文件结构

```
vtune-profile-tool-script/
├── pbgz_wrapper.sh          # 主要的包装脚本
├── README.md                # 本说明文档
└── vtune_results/           # VTune分析结果目录（运行时创建）
```

## 快速开始

### 1. 基本使用

```bash
# 分析当前目录下的 mydata 文件夹
./pbgz_wrapper.sh ./mydata

# 指定输出文件名
./pbgz_wrapper.sh ./mydata my_archive.pbgz

# 传递额外的 pbgz 选项
./pbgz_wrapper.sh ./mydata test.pbgz -t 4 -v
```

### 2. 前置条件

- **pbgz 二进制文件**: 脚本会自动查找 `../release-release-profile/bin/pbgz`
- **Intel VTune Profiler**: 可选，但推荐安装以获得性能分析功能
- **Bash 4.0+**: 脚本使用了一些现代 bash 特性

### 3. 安装 VTune Profiler

```bash
# Ubuntu/Debian
sudo apt update
sudo apt install intel-oneapi-vtune

# 或者从 Intel 官网下载
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/vtune-profiler.html
```

## 详细使用说明

### 脚本参数

```bash
用法: ./pbgz_wrapper.sh <源目录> [输出文件] [pbgz选项]

参数:
    源目录           要压缩的目录路径（必需）
    输出文件         输出文件名（可选，默认: test.pbgz）
    pbgz选项       传递给pbgz的额外选项（可选）

示例:
    ./pbgz_wrapper.sh ./mydata
    ./pbgz_wrapper.sh ./mydata my_archive.pbgz
    ./pbgz_wrapper.sh ./mydata test.pbgz -t 4 -v
```

### 环境变量

| 变量名 | 默认值 | 说明 |
|--------|--------|------|
| `PBGZ_BIN` | `../release-release-profile/bin/pbgz` | pbgz 可执行文件路径 |
| `VTUNE_RESULTS_DIR` | `./vtune_results` | VTune 结果存储目录 |

### VTune 分析类型

脚本会自动运行以下 VTune 分析：

1. **Hotspots Analysis** - 识别 CPU 热点函数
2. **Performance Snapshot** - 系统性能概览
3. **Memory Access Analysis** - 内存访问模式分析
4. **HPC Performance** - 高性能计算特性分析

## 输出结果

### 1. 压缩结果
- 位置: `vtune_results/<输出文件名>`
- 显示原始大小和压缩后大小对比

### 2. VTune 分析结果
- 位置: `vtune_results/vtune_<时间戳>/`
- 包含:
  - 各分析类型的原始数据
  - HTML 格式的报告文件
  - 可用 VTune GUI 打开查看

### 3. 查看结果

```bash
# 使用 VTune GUI 查看热点分析结果
vtune-gui vtune_results/vtune_20231018_160842/hotspots

# 直接在浏览器中查看 HTML 报告
firefox vtune_results/vtune_20231018_160842/hotspots/hotspots_report.html
```

## 使用示例

### 示例 1: 基本性能分析

```bash
# 创建测试数据
mkdir -p test_data
echo "Hello World" > test_data/test.txt
echo "Test data for pbgz" > test_data/data.txt

# 运行分析
./pbgz_wrapper.sh ./test_data
```

### 示例 2: 自定义输出和选项

```bash
# 使用多线程和详细输出
./pbgz_wrapper.sh ./test_data custom_output.pbgz -t 8 -v
```

### 示例 3: 分析大型项目

```bash
# 分析项目源代码
./pbgz_wrapper.sh ./src project_source.pbgz
```

## 故障排除

### 常见问题

#### 1. pbgz 二进制文件未找到
```
[ERROR] pbgz binary not found at: ../release-release-profile/bin/pbgz
```
**解决方案**: 
- 确保 pbgz 已经编译: `./build-release-profile.sh`
- 或设置环境变量: `export PBGZ_BIN=/path/to/pbgz`

#### 2. VTune 命令未找到
```
[WARN] vtune command not found in PATH
```
**解决方案**:
- 安装 Intel VTune Profiler
- 或将 VTune 添加到 PATH: `export PATH=/opt/intel/oneapi/vtune/latest/bin:$PATH`

#### 3. 权限问题
```
[ERROR] Source directory is not readable: /path/to/dir
```
**解决方案**:
- 检查目录权限: `ls -la /path/to/dir`
- 确保有读取权限

#### 4. 磁盘空间不足
```
[ERROR] pbgz compression failed
```
**解决方案**:
- 检查可用磁盘空间: `df -h`
- 清理临时文件: `rm -rf vtune_results/*`

### 调试模式

如果脚本运行出现问题，可以手动调试：

```bash
# 手动运行 pbgz 命令
tar -cf - ./test_data | ../release-release-profile/bin/pbgz -d -c > test.pbgz

# 手动运行 VTune 分析
vtune -collect hotspots -- ./pbgz_wrapper.sh ./test_data
```

## 性能优化建议

基于 VTune 分析结果，通常可以关注以下优化方向：

### 1. CPU 优化
- 查看热点函数，优化算法逻辑
- 检查指令级并行度 (IPC)
- 优化循环结构

### 2. 内存优化
- 减少缓存未命中
- 优化内存访问模式
- 检查内存分配策略

### 3. 并行化优化
- 利用多核处理器
- 减少线程同步开销
- 优化负载均衡

### 4. I/O 优化
- 使用异步 I/O
- 批量处理数据
- 优化文件访问模式

## 高级用法

### 1. 自定义分析脚本

可以修改 `pbgz_wrapper.sh` 来添加自定义的 VTune 分析：

```bash
# 在 run_vtune_profiling 函数中添加自定义分析
local analyses=("hotspots" "performance-snapshot" "memory-access" "hpc-performance" "uarch-exploration")
```

### 2. 批量分析

创建批量分析脚本：

```bash
#!/bin/bash
# batch_analysis.sh
for dir in ./test_data1 ./test_data2 ./test_data3; do
    echo "Analyzing $dir..."
    ./pbgz_wrapper.sh "$dir" "${dir##*/}.pbgz"
done
```

### 3. 持续集成

在 CI/CD 流水线中使用：

```yaml
# .gitlab-ci.yml 示例
performance_test:
  stage: test
  script:
    - ./vtune-profile-tool-script/pbgz_wrapper.sh ./test_data
  artifacts:
    reports:
      junit: vtune_results/vtune_*/hotspots_report.html
```

## 贡献指南

如果您想改进这个工具：

1. Fork 项目
2. 创建功能分支
3. 提交更改
4. 创建 Pull Request

## 许可证

本工具遵循与 pbgz 项目相同的许可证。

## 联系信息

如有问题或建议，请通过以下方式联系：

- 项目 Issues: [GitHub Issues](https://github.com/PBGZip/pbgz/issues)
- 邮件: 项目维护者

---

**注意**: 使用 VTune Profiler 需要有效的 Intel CPU 和适当的系统权限。在某些系统上可能需要管理员权限来安装采样驱动程序。
