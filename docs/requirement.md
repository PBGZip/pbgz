设计思路：
    1. PBGZ 是一种开放、灵活的，面向超大型数据压缩的存储格式，满足：
         - 支持大文件的随机读
         - 支持大文件的写入确定性，且确保每次同一压缩后的数据都完全相同（不会因多线程压缩结果顺序不一致而产生压缩文件不一致）,依赖API 和文件格式设计。提前先分配
         - 支持大文件的增量更新
         - 支持大文件直接 cat。
         - 部分可损毁情况下，解压时，损毁数据可跳过。
         - 允许压缩数据块存在元信息描述，从而支持多种数据块的读写算法。
         - 每个数据块本身都需要校验。
         - 支持流式数据的压缩写入与读出（比如来自stdin, 不清楚数据最总大小）

    2. PBGZ 文件设计：  可损毁部分标记为 （CanDamage）
      
    FILE HEADER:
        magic: "PBGZ" (4 bytes)   （CanDamage）
        version: 1 (1 byte)       （CanDamage）
        origin format: （2 bytes）, fastq? bam？ cram? tar？ generic?  （CanDamage）


    Desc Block:   （CanDamage, optional） if it appears, following  data block's options is defaulted by this desc_block.
        - magic desc_block: 0xSB (1 bytes)   
        - u32 dirblock meta json string size ,4 bytes  （CanDamage,  read policy: search \0 in following bytes）
        - meta json string , contains: from which file? or stdin? or ....       （* CanDamage, if damage , we can manually fix this content, meta json must contains size, Integrity check method, check data,should not use any compression, to make the data more recoverable）

           if follow datablock is fastq {
            - fastq specific meta json, contains: 
            {
                "writer": "pbgz-writer-v1", (required)
                "fixed_len":  xx, (if xx==0 , means variable length) (required)
                "reference_genome_url": "https://xxxx.com/genome.fa", (optional)
                "hash_method": "farmhash", (default: farmhash, optional)
                "reference_genome_hash": "abc123def456...", (optional)
                "cooked_reference_method": "fast", (default: fast, optional)
                "cooked_reference_hash": "xyz789...", (optional if reference_genome_hash is provided)
            }
           } else if follow datablock is bam {
            - bam specific meta json, contains: reference genome, read count, etc.
           } else if follow datablock is tar {
            - tar specific meta json, contains: file list, etc.
           } else if follow datablock is generic {
            - generic specific meta json, contains: data type, etc.
           } else {
            - unknown format
           }

        - \0  ;  （CanDamage ）
        - farmhash checksum of desc_block, 8 bytes  （CanDamage）


    Data Block:  （CanDamage, optional）
        - magic data_block: 0xDB (1 bytes)    
        - content data json meta 
            - meta json string , contains: 
                {
                    "data_type": "reference"
                    "compression_method": "fc", (default: fc, optional)
                    "blocksize": 123456, (size of the data in bytes, required)
                    "reference_genome_url": "https://xxxx.com/genome.fa", (optional)
                    "reference_genome_hash": "abc123def456...", (optional)
                    "cooked_reference_hash": "xyz789...", (required)
                    "cooked_reference_method": "fast", (default: fast, optional)
                }
                or
                {
                    "data_type": "fastq", (default: fastq, optional)
                    "fixed_len": xx, (if xx==0 , means variable length)
                    "compression_method": "fc-bcm-fc-fc", (default: fc-bcm-fc-fc, optional)
                    "blocksize": [ 12, 34, 56, 78], (size of compressed seqid block, read count, etc., required)
                    "reference_genome_url": "https://xxxx.com/genome.fa", (optional)
                    "hash_method": "farmhash", (default: farmhash, optional)
                    "reference_genome_hash": "abc123def456...", (optional)
                } 
                or 
                {
                    "data_type": "paired-fastq", (default: paired-fastq, optional)
                    "fixed_len": xx, (if xx==0 , means variable length)
                    "compression_method": "fc-bcm-fc-fc", (default: fc-bcm-fc-fc, optional)
                    "data_size": [ 12, 34, 56, 78], (size of compressed seqid block, read count, etc., required)
                    "reference_genome_url": "https://xxxx.com/genome.fa", (optional)
                    "hash_method": "farmhash", (default: farmhash, optional)
                    "reference_genome_hash": "abc123def456...", (optional)
                }
        - \0  ;  （CanDamage ）
        - [ u32 * 4 block size , if data type is fastq ]
          [ u32 * 1 block size , if data type is reference]
        - data block content, compressed data, size is blocksize bytes,  （CanDamage）
        - farmhash checksum of data_block, 8 bytes  （CanDamage）
        - origin content checksum of data _block, 8 bytes  （CanDamage）

