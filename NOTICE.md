# PBGZ 开源声明

## 项目许可证

本项目采用 MIT 许可证，详细的许可证内容请参见 [LICENSE](LICENSE) 文件。

## 第三方组件声明

### 项目概述
PBGZ (Parallel BGZF) 是一个并行 BGZF 文件压缩/解压缩库，主要用于生物信息学领域的高性能数据处理。

### 第三方依赖库

本项目使用了以下开源第三方组件，特此致谢：

#### 1. HTSlib (版本 1.12)
- **许可证**: BSD 3-Clause License
- **用途**: 处理 HTS (High-Throughput Sequencing) 数据格式
- **源码**: https://github.com/samtools/htslib
- **修改**: 本项目包含 HTSlib 的部分代码，用于 BGZF 格式支持

#### 2. FarmHash
- **许可证**: MIT License
- **用途**: 哈希函数实现
- **源码**: https://github.com/google/farmhash
- **修改**: 未修改，直接用于文件完整性校验

#### 3. SAIS (Suffix Array Induced Sorting)
- **许可证**: BSD License
- **用途**: 后缀数组排序算法
- **源码**: 包含在 src/coder/bcm_libsais.cpp 中
- **修改**: 集成到项目编码器中

#### 4. ZLIB (版本 1.2.11)
- **许可证**: zlib License
- **用途**: 数据压缩库
- **源码**: https://zlib.net/

#### 5. BZIP2
- **许可证**: BSD-style License
- **用途**: 压缩算法支持
- **源码**: https://sourceware.org/bzip2/

#### 6. LIBDEFLATE (版本 1.7)
- **许可证**: MIT License
- **用途**: 高性能 DEFLATE 压缩/解压缩
- **源码**: https://github.com/ebiggers/libdeflate

#### 7. ZSTD (版本 1.5.0)
- **许可证**: BSD + GPLv2 Dual License
- **用途**: Zstandard 压缩算法支持
- **源码**: https://github.com/facebook/zstd

#### 8. JSONCPP (版本 1.9.4)
- **许可证**: MIT License
- **用途**: JSON 解析和生成
- **源码**: https://github.com/open-source-parsers/jsoncpp

#### 9. Google Test
- **许可证**: BSD 3-Clause License
- **用途**: 单元测试框架
- **源码**: https://github.com/google/googletest

#### 10. SPINLOCK 实现
- **许可证**: various (各实现有不同许可证)
- **用途**: 多种自旋锁实现
- **源码**: 包含在 src/spinlock/ 目录下

### 构建工具
本项目使用了以下构建工具：
- **Autoconf** (版本 2.69) - GPL v3
- **Automake** (版本 1.17) - GPL v3
- **Libtool** (版本 2.4.6) - GPL v3
- **NASM** (版本 2.14.02) - BSD 2-Clause
- **YASM** (版本 1.3.0) - BSD 2-Clause

## 许可证兼容性

本项目选择的 MIT 许可证与所有使用的第三方组件许可证兼容：

- MIT 许可证与 BSD 系列、zlib License、Apache 2.0 等许可证完全兼容
- 与 GPL 系列许可证兼容（本代码不包含 GPL 代码）
- 所有依赖项均为宽松许可证，可以进行商业使用

## 版权声明

### 主要版权
- PBGZ 项目版权所有 (C) 2021
- 详见源文件中的具体版权声明

### 第三方版权
各第三方组件的版权信息请参考其源码包中的 COPYRIGHT 或 LICENSE 文件。

## 免责声明

本软件按"原样"提供，不提供任何明示或暗示的保证。在任何情况下，作者或版权持有人均不对任何索赔、损害或其他责任负责，无论是合同、侵权或其他方面的责任。

## 联系信息

如有关于许可证或第三方使用的疑问，请通过以下方式联系：

- 项目主页: https://github.com/PBGZip/pbgz
- 问题报告: https://github.com/PBGZip/pbgz/issues

## 完整的许可证文本

### MIT 许可证

```
MIT License

Copyright (c) 2021 PBGZ

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

*本文件最后更新时间：2026-01-14*
