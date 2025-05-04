[English](README.en.md)

# PBGZ

- **P**hytium & **B**IG Collaborative Project: A highly competitive **G**enomic and large-scale data **Z**ip tool developed by researchers from Phytium and BIG.  
- **P**etabyte-Scale Capability: Designed for compressing **G**enomic data or generic files at **Z**ip-level efficiency.  

# Key Advantages  

- **Ultra-High Compression Ratio**: Up to **8x better than gzip**.  
- **Blazing-Fast Speed**: Achieve **13x faster compression/decompression** compared to gzip.  
- **Award-Winning Algorithm**: Powered by **FC Algorithm** - ranked **2nd** in the balanced track of the 2021 Global Data Compression Competition.  
- **Cross-Platform Support**: Optimized for x86_64 (Intel/AMD) and ARMv8 (Phytium) CPUs; compatible with Linux (Kylin, Ubuntu) and macOS.  
- **Genomic-First Design**: Specialized for genomic data while supporting generic binary files.  
- **Truly Open Source**: Released under **MIT License** - free for commercial use, modification, and distribution.  

![GDCC 2021 Performance](image/gdcc2021.png)  

## Benchmark Comparison  

| Metric/Feature          | Platform | GZip (.gz) | Enancio (.ora)                      | PBGZ                             |
|-------------------------|----------|------------|--------------------------------------|----------------------------------|
| **Compression Ratio**   | NovaSeq  | 17.01%     | 2.7%                                | **2.4%**                        |
|                         | T7       | 35.51%     | ---                                  | **17.03%**                      |
| **Compression Speed**   | NovaSeq  | 19.54MB/s  | 240MB/s (Whitepaper)                | **265.84MB/s**                  |
|                         | T7       | 12.87MB/s  | ---                                  | **144MB/s**                     |
| **No Dict. Required**   |          | Yes        | No (600MB dict. required)            | Yes (Uses user's FASTA directly)|
| **Multi-Species Support**|          | ---        | Limited by 600MB dict.              | **Full Support**                |
| **Architectures**        |          | x86_64/ARMv8 | x86_64 only                        | x86_64/**ARMv8**               |
| **Ref-Free Decompression**|        | Yes        | No (Requires ref bin)                | Yes                             |
| **MD5 Verification**     |          | No         | Yes                                  | Yes                             |
| **Concatenation Support**|         | Yes        | No                                   | Yes                             |

# Getting Started  

## 1. Installation  

```bash
git clone https://github.com/PBGZip/pbgz.git
cd pbgz
./3rdparty/build.all.gcc.sh
./build-release.sh
```

## 2. Usage  

```bash
Usage: pbgz [OPTION] [FILE]

Mandatory arguments to long options are required for short options too.

-d, --decompress <file.nz> Decompress specified file
-z, --gz Output decompressed file in .gz format
-o, --outfile <file> Specify output filename
-O, --outdir <dir> Specify output directory
-f, --force Overwrite output files
-r, --reference <FASTA> Reference FASTA file (raw .fasta only; .fasta.gz NOT allowed)
-n, --refunpack Embed reference into compressed file (required for decompression)
-t, --threads <N> Number of threads (default: all cores)
-l, --level <1-3> Compression level (1=fast, 2=default, 3=best)
-e, --remove Delete source after successful compression
-h, --help Show help

Compression Example
pbgz human.fq.gz -o human.fq.gz.nz -r ucsc.hg19.fa

Decompression Example
pbgz -d human.fq.gz.nz

```

# Roadmap  

## Pre-v1.0 Tasks  

1. **Critical Fixes**  
   a. Special handling for **first-gen sequencing data** (variable-length reads)  
   b. Cross-genome compression stability testing  

2. **Core Integrity Features**  
   a. Checksum validation for ref bin data  
   b. Store FASTA MD5 in metadata to prevent wrong-reference decompression  
   c. Metadata redundancy across file regions  

3. **Format Enhancements**  
   a. Support concatenated .pbgz files  
   b. Optimized paired-end FASTQ handling  

4. **Future-Proof Design**  
   - Ensure format extensibility for new features  

5. **Scalability Tests**  
   - Ultra-large genome testing (e.g., wheat genome)  

6. **Pipeline Integration**  
   - Support stdin/stdout streams and SSH-based workflows  

7. **Code Compliance**  
   a. Restore license headers  
   b. Add English comments (keep Chinese where clearer)  

8. **Documentation**  
   - Complete build guides and technical specs  

## Post-v1.0 Plans  

1. **Extended Format Support**  
   - Third-gen sequencing data (PacBio/Nanopore)  
   - HDF5 I/O integration  

2. **Domain Expansion**  
   - CT (Computed Tomography) data support  

3. **Library Development**  
   - Create libpbgz.so with Python/Java/C++ APIs  
