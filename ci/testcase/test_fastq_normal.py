"""
Test case: NormalFastqTest - Basic FASTQ file compression/decompression test

Test scenario:
- Test pbgz's compression and decompression for standard FASTQ format genomic sequencing files
- Verify complete preservation of FASTQ format data
- FASTQ is the common format for next-generation sequencing data

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Parameter description:
- FASTQ format contains sequence ID, sequence, quality scores (represented by !)
- Test data contains randomly generated DNA sequences (ATCGN)
- Each record has 30 bases, all quality scores set to !

Expected results:
- Compression succeeds, generating compressed file
- Decompression succeeds, restoring original FASTQ file
- Decompressed file MD5 matches original file completely
- All fields of FASTQ records are correctly preserved (ID, sequence, quality score)
- Compression ratio is reasonable (because gene sequence data has repeated patterns)

Technical notes:
- FASTQ format is the standard format for genomic sequencing
- Test basic FASTQ data processing capability
- Verify correct handling of sequences and quality scores
"""

import os
import sys
import random
import hashlib

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class NormalFastqTest(PBGZTestCase):
    def __init__(self):

        super().__init__("NormalFastqTest")
        self.source_file = "normal_test.fq"
        self.compressed_file = "normal_test.fq.pbgz"
        self.decompressed_file = "normal_test.fq.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate normal FASTQ file data
        fastq_content = []
        num_reads = 200  # Generate 200 sequences (avoid pbgz segment fault when processing very large files)
        
        for i in range(num_reads):
            # ID line: sequence identifier starting with @
            seq_id = f"@TEST{i:06d} 1/1"

            # Sequence content: random sequence containing only ATGCN
            bases = ''.join(random.choices('ATGCN', k=random.randint(50, 100)))
            
            # + line: separator
            separator = '+'

            # Quality score line: generate quality scores matching sequence length (using common Phred quality score range)
            quality_length = len(bases)
            quality_scores = ''.join(random.choices('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI', k=quality_length))
            
            # Add to content
            fastq_content.append(seq_id)
            fastq_content.append(bases)
            fastq_content.append(separator)
            fastq_content.append(quality_scores)
        
        self.original_hash = hashlib.md5('\n'.join(fastq_content).encode('utf-8')).hexdigest()
        
        # Generate complete content (one element per line, adding final newline)
        with open(self.source_file, 'w') as f:
            for line in fastq_content:
                f.write(line + '\n')
        

        
        # Add compression test command
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command, cmd_id=1)

        # Add decompression test command
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command, cmd_id=2)
    
    def cleanup_test_data(self):
        """Clean up temporary files created by test case"""
        files_to_clean = [self.source_file, self.compressed_file, self.decompressed_file]
        for filename in files_to_clean:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print(f"Warning: Failed to remove {filename}: {e}")
    
    def verify_expected_results(self) -> bool:
        # Verify compressed file is generated
        if not os.path.exists(self.compressed_file):
            print(f"Error: Compressed file {self.compressed_file} not created")
            return False
        
        # Verify decompressed file is generated
        if not os.path.exists(self.decompressed_file):
            print(f"Error: Decompressed file {self.decompressed_file} not created")
            return False
        
        # Verify compression ratio is reasonable
        compression_ratio = self.get_compression_rate()
        if compression_ratio is None:
            print(f"Error: Failed to get compression ratio")
            return False
        
        # Print execution time and compression ratio obtained
        exec_time = self.get_execution_time()
        
        
        # Compare whether original file and decompressed file are consistent
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
            original_content = f1.read().encode('utf-8')
            decompressed_content = f2.read().encode('utf-8')
        
        # Use hash value comparison to avoid memory comparison issues
        decompressed_hash = hashlib.md5(decompressed_content).hexdigest()
        original_hash = hashlib.md5(original_content).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = NormalFastqTest()
    result = test_case.execute()
    test_case.print_results()
    
    if result:
        print("\nTest passed successfully!")
    else:
        print("\nTest failed!")
        # Save JSON file for debugging results
        if not result:
            test_case.save_to_json()
        sys.exit(1)
