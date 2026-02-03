# Changelog

This document records important changes to the PBGZ project.

## [2.1.0] - 2026-01-26

### Major Features
- Support for secondary command functionality

### Important Fixes
- Fixed decompression failure issue in reference genome compression without packaging scenario
- Fixed abnormal exit issue after compression completion

## [2.0.0] - 2025-12-26

### Major Features
- Support for reference genome compression and decompression functionality
- Support for pipeline reading and writing functionality

### Architecture Refactoring
- Redesigned pbgz file format
- Optimized project structure, adjusted code to src directory

## [1.0.0] - 2025-09-26

### Major Features
- Initial version release
- Core compression engine based on FC high-efficiency compression algorithm
- Support for x86_64 and ARMv8 architectures
- Support for Linux and Mac OSX operating systems

### Core Advantages
- Extreme compression ratio: up to 8x that of gzip
- Extreme compression/decompression speed: up to 13x that of gzip
- Support for multi-species reference genomes
- Decompression does not require reference genome

---

*This document records important features and issue fixes of the project. For detailed changes, please refer to the git commit history.*
