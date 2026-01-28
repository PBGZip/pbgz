# PBGZ 

- 可以理解为： **P**hytium & **B**IG 的研究人员一起发起的一个极具竞争力的高速.高压缩比大型 **G**enomic数据或大型通用数据的压缩(**Z**ip)工具。
- 也可以理解为：支持PB级别的 **G**enomic数据或大型通用数据的压缩(**Z**ip)工具。

# PBGZ 的优势

- 极致压缩率：最高可达 gzip 的8倍。
- 极致加解压速度：最高可达 gzip 的13倍。
- 先进算法：采用 Phytium ResearchLab 在2021 全球压缩竞赛中 平衡组第2名的 FC 高效压缩算法
- 多平台支持：在Intel/AMD等x86_64 和 飞腾等ARMv8体系结构 CPU 上均可高效运行；支持 麒麟.Ubuntu等Linux 和 MAC OSX 操作系统。
- 多种数据支持：针对基因数据专门优化，但同时支持任意二进制文件的压缩；
- 完全自由开放：基于 MIT License 开源，允许任何人.任何组织使用.修改和分发。

![alt text](images/gdcc2021.png)


| 指标/特性              | 平台    | GZip (.gz)    | Enancio (.ora)                               | Ours                             |
|----------------------|---------|---------------|-----------------------------------------------|---------------------------------|
| **压缩率**            | NovaSeq | 17.01%        | 2.7%                                          | 2.4%                            |
|                      | T7      | 35.51%        | ---                                           | 17.03%                          |
| **压缩速度**           | NovaSeq| 19.54MB/s     | 240MB/s (白皮书)                               | **265.84MB/s**                   |
|                       | T7     | 12.87MB/s     | ---                                           | **144MB/s**                     |
| **发布包不需要字典数据** |        | Yes           | No, 发布包必须配600MB的字典                      | Yes (直接使用用户的fasta)           |
| **支持多物种参考基因组** |        | ---           | No, 目前受发布包600MB字典数据影响，实际很难支持多种物种基因数据。 | Yes, 可以广泛支持所有物种  |
| **支持的体系结构**      |        | x86_64，ARMv8  | x86_64                                        | x86_64，ARMv8                    |
| **解压不需要参考基因组** |        | Yes           | No，必须依赖发布包的bin数据                       | Yes                              |
| **支持MD5校验**        |        | No            | Yes                                           | Yes                              |
| **支持多个压缩包直接合并**|       | Yes             | No                                           | Yes                             |

## 安装说明

### 系统要求

- **操作系统**: Linux (Ubuntu, 麒麟等), macOS
- **架构**: x86_64 (Intel/AMD), ARMv8 (飞腾等)
- **编译器**: GCC 7.0+ 或 Clang 8.0+
- **CMake**: 版本 3.10 或更高
- **内存**: 建议 4GB 以上内存用于编译

### 依赖库

PBGZ 需要以下第三方库支持：

- **Intel ISA-L**: 用于优化压缩算法
- **libdeflate**: 高性能压缩库
- **HTSlib**: 基因数据处理库
- **JSONCPP**: JSON 处理库
- **Zstd**: Zstandard 压缩库
- **Bzip2**: bzip2 压缩库
- **Zlib**: zlib 压缩库

### 快速安装

#### 1. 克隆仓库

```bash
git clone https://github.com/PBGZip/pbgz.git
cd pbgz
```

#### 2. 构建第三方依赖库

```bash
# 首先构建所有第三方依赖库
cd 3rd_party
./build.all.gcc.sh
cd ..
```

#### 3. 编译主程序

```bash
# 编译 Release 版本（推荐用于生产环境）
./build-release.sh

# 或者编译 Debug 版本（用于开发调试）
./build-debug.sh
```

#### 4. 验证安装

```bash
# 检查依赖库
./scripts/check_dependencies.sh

# 运行测试
./test/pbgz_test
```

### 预编译版本

如果您不想从源码编译，可以使用预编译版本：

```bash
# 下载预编译版本
wget https://github.com/PBGZip/pbgz/releases/latest/download/pbgz-linux-x86_64.tar.gz

# 解压
tar -xzf pbgz-linux-x86_64.tar.gz

# 添加到 PATH
export PATH=$PWD/pbgz/bin:$PATH
```

## 使用指南

### 基本用法

PBGZ 采用子命令设计，支持三种主要操作：压缩、解压和创建索引。

```bash
# 查看帮助
pbgz help                    # 显示通用帮助
pbgz help compress           # 显示压缩命令帮助
pbgz help decompress         # 显示解压命令帮助

# 查看版本
pbgz version
```

#### 压缩文件

```bash
# 基本压缩
pbgz compress input.fastq -o output.pbgz

# 使用参考基因组（推荐，可获得更好的压缩比）
pbgz compress input.fastq -o output.pbgz -r reference.fa

# 指定输出目录
pbgz compress input.fastq -O /output/directory/

# 指定线程数和压缩级别
pbgz compress input.fastq -o output.pbgz -t 8 -l 6
```

#### 解压文件

```bash
# 基本解压
pbgz decompress input.pbgz

# 指定输出文件名
pbgz decompress input.pbgz -o output.fastq

# 解压为gzip格式
pbgz decompress input.pbgz -z

# 使用参考基因组（如果需要）
pbgz decompress input.pbgz -r reference.fa
```

### 命令行参数说明

PBGZ 使用子命令格式，主要参数如下：

#### 压缩参数 (compress)
- `-o, --outfile`: 指定输出文件名
- `-O, --outdir`: 指定输出目录
- `-r, --reference`: 指定参考基因组文件（FASTA格式，基因数据压缩推荐）
- `-t, --threads`: 指定线程数（默认为CPU核心数）
- `-l, --level`: 压缩级别（1-9，1最快，9压缩比最好，默认5）
- `-e, --remove`: 压缩成功后删除原文件
- `-f, --force`: 强制覆盖输出文件
- `-i, --index`: 压缩时创建索引文件

#### 解压参数 (decompress)
- `-o, --outfile`: 指定输出文件名
- `-O, --outdir`: 指定输出目录
- `-r, --reference`: 指定参考基因组文件（某些压缩文件需要）
- `-p, --position`: 指定基因位置范围（如chr1:1000-2000）
- `-z, --gz`: 解压为gzip格式
- `-f, --force`: 强制覆盖输出文件

### 使用示例

#### 示例1：压缩 FASTQ 文件


```bash
# 压缩基因测序数据
pbgz compress sample_NovaSeq.fastq -o sample_NovaSeq.pbgz -r hg38.fa -l 9 -t 16
```

#### 示例2：解压缩

```bash
# 解压并验证
pbgz decompress sample_NovaSeq.pbgz -o sample_NovaSeq.fastq -t 16
```

### 性能优化建议

#### 压缩优化

1. **使用参考基因组**: 对于基因数据，使用参考基因组可以显著提高压缩率
2. **选择合适的压缩级别**: 
   - 级别 1-3: 速度快，压缩率一般
   - 级别 4-6: 平衡速度和压缩率（推荐）
   - 级别 7-9: 压缩率最高，速度较慢
3. **多线程**: 根据CPU核心数设置合适的线程数
4. **块大小**: 对于大文件，可以增加块大小以提高压缩效率

#### 解压优化

1. **多线程解压**: 使用多线程可以显著提高解压速度
2. **验证选项**: 在生产环境中建议开启 MD5 验证
3. **内存使用**: 解压过程会使用一定内存，建议系统内存充足

### 故障排除

#### 常见问题

1. **依赖库缺失**
   ```bash
   # 检查依赖
   ./scripts/check_dependencies.sh
   
   # 重新构建依赖
   cd 3rd_party && ./build.intel.sh
   ```

2. **运行时错误**
   ```bash
   # 检查库文件路径
   export LD_LIBRARY_PATH=/path/to/pbgz/3rd_party/release/lib:$LD_LIBRARY_PATH
   
   # 使用日志模式查看错误信息
   pbgz compress input.txt -o output.pbgz -g 2
   ```


#### 调试模式

启用详细输出进行问题诊断：

```bash
pbgz compress input.fastq -o output.pbgz -g 1
```

### 高级功能

#### 流式压缩

PBGZ 支持流式压缩，适合处理大文件或管道操作：

```bash
# 从管道压缩
cat large_file.fastq | pbgz compress -o compressed.pbgz

# 解压到管道
pbgz decompress compressed.pbgz -o - | grep "sequence"
```

## 贡献指南

我们欢迎社区贡献！请查看 [CONTRIBUTING.md](CONTRIBUTING.md) 了解详细信息。

## 许可证

本项目基于 MIT 许可证开源。详见 [LICENSE](LICENSE) 文件。

## 致谢

- Intel ISA-L 库用于压缩算法优化
- HTSlib 库用于基因数据处理
- FC 算法在 2021 全球压缩竞赛中获得平衡组第2名

## 联系我们

- 项目主页: https://github.com/PBGZip/pbgz
- 问题报告: https://github.com/PBGZip/pbgz/issues
- 邮箱: wenjinyang2729@phytium.com.cn
