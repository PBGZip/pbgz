"""
Test case: FastqExtraSeparatorTest - FASTQ file extra separator handling test

Test scenario:
- Test pbgz handling FASTQ files containing extra separators
- Verify pbgz's fault tolerance for non-standard FASTQ formats
- Test separator processing data integrity

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Parameter description:
- FASTQ format requires strict @ and + separators
- Test extreme cases containing extra separators
- Verify file format fault tolerance

Expected results:
- Compression succeeds, handling non-standard format
- Decompression succeeds, data completely preserved
- MD5 verification matches completely
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqExtraSeparatorTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqExtraSeparatorTest")
        self.source_file = "extra_separator.fq"
        self.compressed_file = "extra_separator.fq.pbgz"
        self.decompressed_file = "extra_separator.fq.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate FASTQ file with more separator lines than ID lines
        fastq_content = []
        num_reads = 200  # Generate 200 sequences

        for i in range(num_reads):
            # ID line: sequence identifier starting with @
            seq_id = f"@TEST{i:06d} 1/1"

            # Sequence content: random sequence containing only ATGCN
            bases = ''.join(random.choices('ATGCN', k=random.randint(50, 100)))
            
            # + line: separator - deliberately create separator count anomaly at record 50 (extra separator)
            if i == 50:
                # Separator anomaly: extra separator line
                separator = '+'
                fastq_content.extend([seq_id, bases, separator])
                # Extra separator line inserted here
                fastq_content.append('+extra')
                # Add quality score line
                quality_length = len(bases)
                quality = ''.join(random.choice('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI')
                                 for _ in range(quality_length))
                fastq_content.append(quality)
            else:
                # Normal + separator
                separator = '+'
                fastq_content.extend([seq_id, bases, separator])
                quality_length = len(bases)
                quality = ''.join(random.choice('\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                                 for _ in range(quality_length))
                fastq_content.append(quality)
        
        full_content = 'n'.join(fastq_content) + 'n'  # Normal FASTQ file ends with newline
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

        # Print execution time and compression ratio
        exec_time = self.get_execution_time()


        # Verify compression ratio is valid
        if compression_ratio <= 0:
            print(f"Error: Invalid compression ratio: {compression_ratio:.2f}%")
            return False

        # Compare original file and decompressed file for consistency
        with open(self.source_file, 'r', encoding='utf-8') as f1, open(self.decompressed_file, 'r', encoding='utf-8') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()

        # Use hash value comparison
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = FastqExtraSeparatorTest()
    result = test_case.execute()
    test_case.print_results()
    
    if result:
        print("\nTest passed successfully!")
    else:
        print("\\nTest failed!")
        # Save JSON file for debugging results
        if not result:
            test_case.save_to_json()
        sys.exit(1)
