"""
Test case: FastqInvalidBaseTest - FASTQ file invalid base handling test

Test scenario:
- Test FASTQ files containing invalid bases (non-ATCGN characters)
- Verify pbgz's fault tolerance capability and handling strategy for invalid bases
- Test data format validation and exception handling

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Parameter description:
- FASTQ standard bases are A, T, C, G, N, test illegal characters like Z
- Verify pbgz's handling of abnormal data
- Test balance between format strictness and fault tolerance

Expected results:
- Compression succeeds, handling invalid base data
- Decompression succeeds, preserving original data (including illegal characters)
- MD5 verification matches completely
"""

import os
import sys
import random
import hashlib

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqInvalidBaseTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqInvalidBaseTest")
        self.source_file = "invalid_base_test.fq"
        self.compressed_file = "invalid_base_test.fq.pbgz"
        self.decompressed_file = "invalid_base_test.fq.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate FASTQ file containing invalid bases
        fastq_content = []
        num_reads = 50

        for i in range(num_reads):
            # ID line
            seq_id = f"@INVALID{i:06d} 1/1"

            # Sequence content: insert some invalid characters into random sequence
            valid_bases = ''.join(random.choices('ATGCN', k=45))
            # Insert invalid characters into sequence (non-ATGCN)
            if i % 3 == 0:
                seq = valid_bases + 'XZ'  # Insert X and Z
            elif i % 3 == 1:
                seq = valid_bases + 'KLM'  # Insert KLM
            else:
                seq = valid_bases + 'OPQR'  # Insert OPQR

            # + line
            separator = '+'

            # Quality score line
            quality_length = len(seq)
            quality = ''.join(random.choice('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                             for _ in range(quality_length))
            
            fastq_content.extend([seq_id, seq, separator, quality])
        
        full_content = '\n'.join(fastq_content) + '\n'
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        
        
        # Add compression test command
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command)

        # Add decompression test command
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command)
    
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
        
        
        # Compare original file and decompressed file for consistency
        with open(self.source_file, 'r', encoding='utf-8') as f1, open(self.decompressed_file, 'r', encoding='utf-8') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = FastqInvalidBaseTest()
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