# PBGZ Open Source Notice

## Project License

This project is licensed under the MIT License. For detailed license terms, please see the [LICENSE](LICENSE) file.

## Third-Party Components Notice

### Project Overview

- **PBGZ** can be understood as: A highly competitive high-speed, high-compression-ratio large-scale **G**enomic data or large-scale general data compression (**Z**ip) tool initiated by researchers from **Phytium** & **BIG**.
- **PBGZ** can also be understood as: A compression (**Zip**) tool that supports petabyte-level **G**enomic data or large-scale general data.

### Advantages of PBGZ 🚀

- **📦 Extreme Compression Ratio**: Up to 8x higher than gzip
- **⚡ Extreme Compression/Decompression Speed**: Up to 13x faster than gzip
- **🧠 Advanced Algorithm**: Uses the FC high-efficiency compression algorithm that won 2nd place in the balanced category at the 2021 Global Compression Competition by Phytium ResearchLab
- **🌐 Multi-platform Support**: Efficiently runs on both Intel/AMD x86_64 and Phytium ARMv8 architecture CPUs; supports Kylin, Ubuntu, and other Linux as well as MAC OSX operating systems
- **🧬 Multiple Data Support**: Specifically optimized for genomic data while supporting compression of arbitrary binary files
- **🔓 Fully Free and Open**: Based on MIT License open source, allowing any individual or organization to use, modify, and distribute

### Third-Party Dependency Libraries 📚

This project uses the following open-source third-party components, and we hereby express our gratitude 🙏:

#### 1. 🧬 HTSlib (Version 1.12)
- **📜 License**: BSD 3-Clause License
- **🎯 Purpose**: Handling HTS (High-Throughput Sequencing) data formats
- **🔗 Source**: https://github.com/samtools/htslib
- **🔧 Modifications**: This project includes partial code from HTSlib for BGZF format support

#### 2. 🔑 FarmHash
- **📜 License**: MIT License
- **🎯 Purpose**: Hash function implementation
- **🔗 Source**: https://github.com/google/farmhash
- **🔧 Modifications**: Unmodified, directly used for file integrity verification

#### 3. 📊 SAIS (Suffix Array Induced Sorting)
- **📜 License**: BSD License
- **🎯 Purpose**: Suffix array sorting algorithm
- **📂 Source**: Included in src/coder/bcm_libsais.cpp
- **🔧 Modifications**: Integrated into the project encoder

#### 4. 🗜️ ZLIB (Version 1.2.11)
- **📜 License**: zlib License
- **🎯 Purpose**: Data compression library
- **🔗 Source**: https://zlib.net/

#### 5. 📦 BZIP2
- **📜 License**: BSD-style License
- **🎯 Purpose**: Compression algorithm support
- **🔗 Source**: https://sourceware.org/bzip2/

#### 6. ⚡ LIBDEFLATE (Version 1.7)
- **📜 License**: MIT License
- **🎯 Purpose**: High-performance DEFLATE compression/decompression
- **🔗 Source**: https://github.com/ebiggers/libdeflate

#### 7. 🚀 ZSTD (Version 1.5.0)
- **📜 License**: BSD + GPLv2 Dual License
- **🎯 Purpose**: Zstandard compression algorithm support
- **🔗 Source**: https://github.com/facebook/zstd

#### 8. 📝 JSONCPP (Version 1.9.4)
- **📜 License**: MIT License
- **🎯 Purpose**: JSON parsing and generation
- **🔗 Source**: https://github.com/open-source-parsers/jsoncpp

#### 9. 🧪 Google Test
- **📜 License**: BSD 3-Clause License
- **🎯 Purpose**: Unit testing framework
- **🔗 Source**: https://github.com/google/googletest

#### 10. 🔄 SPINLOCK Implementations
- **📜 License**: various (different implementations have different licenses)
- **🎯 Purpose**: Various spinlock implementations
- **📂 Source**: Included in src/spinlock/ directory

### 🔨 Build Tools
This project uses the following build tools:
- **🔧 Autoconf** (Version 2.69) - GPL v3
- **⚙️ Automake** (Version 1.17) - GPL v3
- **🔗 Libtool** (Version 2.4.6) - GPL v3
- **💻 NASM** (Version 2.14.02) - BSD 2-Clause
- **⚡ YASM** (Version 1.3.0) - BSD 2-Clause

## 🔐 License Compatibility

The MIT License chosen for this project is compatible with all used third-party component licenses ✅:

- **✅ MIT License** is fully compatible with BSD series, zlib License, Apache 2.0, and other licenses
- **✅ Compatible** with GPL series licenses (this code does not contain GPL code)
- **🏢 All dependencies** are permissive licenses and can be used commercially

## ©️ Copyright Notice

### 📋 Main Copyright
- **PBGZ Project** Copyright (C) 2025
- See specific copyright notices in source files for details

### 📚 Third-Party Copyright
For copyright information of third-party components, please refer to the COPYRIGHT or LICENSE files in their source packages.

## ⚠️ Disclaimer

This software is provided "as is", without any express or implied warranties. In no event shall the authors or copyright holders be liable for any claim, damages, or other liability, whether in contract, tort, or otherwise, arising from, out of, or in connection with the software or the use or other dealings in the software.

## 📞 Contact Information

For questions regarding licenses or third-party usage, please contact us through the following methods:

- 🏠 **Project homepage**: https://github.com/PBGZip/pbgz
- 🐛 **Issue reporting**: https://github.com/PBGZip/pbgz/issues
- 💬 **Discussions**: https://github.com/PBGZip/pbgz/discussions

## Complete License Text

### MIT License

```
MIT License

Copyright (c) 2025 PBGZip

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

### 📊 Dependency Summary Table

| Dependency | Version | License | Purpose | Status |
|------------|---------|---------|---------|--------|
| HTSlib | 1.12 | BSD 3-Clause | HTS data formats | ✅ Integrated |
| FarmHash | Latest | MIT | Hash functions | ✅ Used |
| SAIS | Embedded | BSD | Sorting algorithm | ✅ Integrated |
| ZLIB | 1.2.11 | zlib | Compression | ✅ Supported |
| BZIP2 | Latest | BSD | Compression | ✅ Supported |
| LIBDEFLATE | 1.7 | MIT | DEFLATE | ✅ Supported |
| ZSTD | 1.5.0 | BSD + GPLv2 | Zstandard | ✅ Supported |
| JSONCPP | 1.9.4 | MIT | JSON parsing | ✅ Used |
| Google Test | Latest | BSD 3-Clause | Testing | ✅ Dev-only |
| SPINLOCK | Multiple | Various | Concurrency | ✅ Integrated |

### 🔍 License Compatibility Matrix

| Our License | BSD | MIT | zlib | Apache 2.0 | GPL | LGPL |
|-------------|-----|-----|------|------------|-----|------|
| **MIT** | ✅ Compatible | ✅ Compatible | ✅ Compatible | ✅ Compatible | ✅ Compatible | ✅ Compatible |

---

*This file was last updated on: January 2026*

---

**📢 Note**: This notice is provided for transparency and compliance purposes. All third-party components are used according to their respective licenses. If you believe any attribution is missing or incorrect, please file an issue on our GitHub repository. 🙏
