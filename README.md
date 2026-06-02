# PBGZ

[English](README.md) | [中文](README.zh-CN.md)

- Can be understood as: A highly competitive, high-speed, high-compression-ratio compression (Zip) tool for large-scale Genomic data or large general-purpose data, jointly launched by researchers from Phytium and BIG.
- Can also be understood as: A compression (Zip) tool supporting PB-level Genomic data or large general-purpose data.

# Advantages of PBGZ

- **Extreme compression ratio**: Up to 8 times that of gzip
- **Extreme compression/decompression speed**: Up to 13 times that of gzip  
- **Advanced algorithm**: Uses the FC high-efficiency compression algorithm from Phytium ResearchLab, which ranked 2nd in the balanced category of the 2021 Global Compression Championship
- **Multi-platform support**: Efficiently runs on Intel/AMD x86_64 and Phytium ARMv8 architecture CPUs; supports Linux (Kylin, Ubuntu, etc.) and macOS operating systems
- **Multiple data support**: Specifically optimized for genomic data, but also supports compression of any binary files
- **Completely free and open**: Open source under MIT License, allowing anyone and any organization to use, modify, and distribute

![alt text](images/gdcc2021.png)


| Metric/Feature                          | Platform  | GZip (.gz)    | Enancio (.ora)                               | Ours                             |
|----------------------------------------|----------|---------------|-----------------------------------------------|---------------------------------|
| **Compression Ratio**                 | NovaSeq   | 17.01%        | 2.7%                                          | 2.4%                            |
|                                        | T7        | 35.51%        | ---                                           | 17.03%                          |
| **Compression Speed**                 | NovaSeq   | 19.54MB/s     | 240MB/s (Whitepaper)                          | **265.84MB/s**                   |
|                                        | T7        | 12.87MB/s     | ---                                           | **144MB/s**                     |
| **Release package doesn't need dictionary data** |           | Yes           | No, release package must be paired with 600MB dictionary | Yes (directly uses user's fasta) |
| **Supports multiple species reference genomes** |           | ---           | No, currently limited by 600MB dictionary data in release package, difficult to support multiple species genomic data | Yes, can widely support all species |
| **Supported architectures**           |           | x86_64, ARMv8 | x86_64                                        | x86_64, ARMv8                    |
| **Decompression doesn't need reference genome** |           | Yes           | No, must rely on release package's bin data   | Yes                              |
| **Supports MD5 verification**         |           | No            | Yes                                           | Yes                              |
| **Supports direct merging of multiple compressed packages** |         | Yes           | No                                            | Yes                             |

## Installation Instructions

### System Requirements

- **Operating System**: Linux (Ubuntu, Kylin, etc.), macOS
- **Architecture**: x86_64 (Intel/AMD), ARMv8 (Phytium, etc.)
- **Compiler**: GCC 7.0+ or Clang 8.0+
- **CMake**: Version 3.10 or higher
- **Memory**: Recommended 4GB+ memory for compilation

### Dependency Libraries

PBGZ requires the following third-party libraries:

- **Intel ISA-L**: For optimizing compression algorithms
- **libdeflate**: High-performance compression library
- **HTSlib**: Genomic data processing library
- **JSONCPP**: JSON processing library
- **Zstd**: Zstandard compression library
- **Bzip2**: bzip2 compression library
- **Zlib**: zlib compression library

### Quick Installation

#### 1. Clone Repository

```bash
git clone https://github.com/PBGZip/pbgz.git
cd pbgz
```

#### 2. Build Third-Party Dependencies

```bash
# First build all third-party dependencies
cd 3rd_party
./build.all.gcc.sh
cd ..
```

#### 3. Compile Main Program

```bash
# Compile Release version (recommended for production environment)
./build-release.sh

# Or compile Debug version (for development and debugging)
./build-debug.sh
```

#### 4. Verify Installation

```bash
# Check dependencies
./scripts/check_dependencies.sh

# Run tests
./test/pbgz_test
```

### Precompiled Version

If you don't want to compile from source, you can use a precompiled version:

```bash
# Download precompiled version
wget https://github.com/PBGZip/pbgz/releases/latest/download/pbgz-linux-x86_64.tar.gz

# Extract
tar -xzf pbgz-linux-x86_64.tar.gz

# Add to PATH
export PATH=$PWD/pbgz/bin:$PATH
```

## Usage Guide

### Basic Usage

PBGZ uses a subcommand design, supporting three main operations: compression, decompression, and index creation.

```bash
# View help
pbgz help                    # Display general help
pbgz help compress           # Display compression command help
pbgz help decompress         # Display decompression command help

# View version
pbgz version
```

#### Compress Files

```bash
# Basic compression
pbgz compress input.fastq -o output.pbgz

# Use reference genome (recommended for better compression ratio)
pbgz compress input.fastq -o output.pbgz -r reference.fa

# Specify output directory
pbgz compress input.fastq -O /output/directory/

# Specify thread count and compression level
pbgz compress input.fastq -o output.pbgz -t 8 -l 6
```

#### Decompress Files

```bash
# Basic decompression
pbgz decompress input.pbgz

# Specify output filename
pbgz decompress input.pbgz -o output.fastq

# Decompress to gzip format
pbgz decompress input.pbgz -z

# Use reference genome (if needed)
pbgz decompress input.pbgz -r reference.fa
```

### Command Line Parameters

PBGZ uses subcommand format with the following main parameters:

#### Compression Parameters (compress)
- `-o, --outfile`: Specify output filename
- `-O, --outdir`: Specify output directory
- `-r, --reference`: Specify reference genome file (FASTA format, recommended for genomic data compression)
- `-t, --threads`: Specify thread count (default is number of CPU cores)
- `-l, --level`: Compression level (1-9, 1 is fastest, 9 is best compression ratio, default 5)
- `-e, --remove`: Delete original file after successful compression
- `-f, --force`: Force overwrite output file
- `-i, --index`: Create index file during compression

#### Decompression Parameters (decompress)
- `-o, --outfile`: Specify output filename
- `-O, --outdir`: Specify output directory
- `-r, --reference`: Specify reference genome file (some compressed files require this)
- `-p, --position`: Specify genomic position range (e.g., chr1:1000-2000)
- `-z, --gz`: Decompress to gzip format
- `-f, --force`: Force overwrite output file

### Usage Examples

#### Example 1: Compress FASTQ Files


```bash
# Compress genomic sequencing data
pbgz compress sample_NovaSeq.fastq -o sample_NovaSeq.pbgz -r hg38.fa -l 9 -t 16
```

#### Example 2: Decompression

```bash
# Decompress and verify
pbgz decompress sample_NovaSeq.pbgz -o sample_NovaSeq.fastq -t 16
```

### Performance Optimization Tips

#### Compression Optimization

1. **Use reference genome**: For genomic data, using a reference genome can significantly improve compression ratio
2. **Choose appropriate compression level**: 
   - Level 1-3: Fast speed, general compression ratio
   - Level 4-6: Balance speed and compression ratio (recommended)
   - Level 7-9: Highest compression ratio, slower speed
3. **Multi-threading**: Set appropriate thread count based on number of CPU cores
4. **Block size**: For large files, you can increase block size to improve compression efficiency

#### Decompression Optimization

1. **Multi-threaded decompression**: Using multiple threads can significantly improve decompression speed
2. **Verification options**: MD5 verification is recommended in production environments
3. **Memory usage**: The decompression process will use some memory, so sufficient system memory is recommended

### Troubleshooting

#### Common Issues

1. **Missing dependencies**
    ```bash
    # Check dependencies
    ./scripts/check_dependencies.sh
    
    # Rebuild dependencies
    cd 3rd_party && ./build.intel.sh
    ```

2. **Runtime errors**
    ```bash
    # Check library file path
    export LD_LIBRARY_PATH=/path/to/pbgz/3rd_party/release/lib:$LD_LIBRARY_PATH
    
    # Use log mode to view error messages
    pbgz compress input.txt -o output.pbgz -g 2
    ```


#### Debug Mode

Enable verbose output for problem diagnosis:

```bash
pbgz compress input.fastq -o output.pbgz -g 1
```

### Advanced Features

#### Streaming Compression

PBGZ supports streaming compression, suitable for processing large files or pipeline operations:

```bash
# Compress from pipeline
cat large_file.fastq | pbgz compress -o compressed.pbgz

# Decompress to pipeline
pbgz decompress compressed.pbgz -o - | grep "sequence"
```

## Contributing Guide

We welcome community contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## License

This project is open source under the MIT License. See [LICENSE](LICENSE) file for details.

## Acknowledgments

- Intel ISA-L library for compression algorithm optimization
- HTSlib library for genomic data processing
- FC algorithm, which ranked 2nd in the balanced category of the 2021 Global Compression Championship

## Contact Us

- Project homepage: https://github.com/PBGZip/pbgz
- Issue reporting: https://github.com/PBGZip/pbgz/issues
- Email: wenjinyang2729@phytium.com.cn