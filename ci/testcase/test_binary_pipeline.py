"""
Test case: BinaryPipelineTest - Binary file multi-stage pipeline processing test

Test scenario:
- Test pbgz's multi-stage processing capability in pipelines
- Verify intermediate processing steps in pipelines
- Test stability in complex pipeline scenarios

Commands used:
1. cat {source_file} | ./release-release/bin/pbgz compress | cat | ./release-release/bin/pbgz decompress -o {decompressed_file}

Expected results:
- Multi-stage pipeline processing succeeds
- Compression and decompression work normally in pipelines
- Intermediate cat command integration works normally
- Final decompressed file matches original file

Technical notes:
- Test complex pipeline combinations
- Verify pbgz stability in pipeline chains
- Test error handling and stream management
"""

import os
import sys
import random
import hashlib

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class BinaryPipelineTest(PBGZTestCase):
    def __init__(self):

        super().__init__("BinaryPipelineTest")
        self.source_file = "pipeline_test_data.bin"
        self.compressed_file = "pipeline_test_data.bin.pbgz"
        self.decompressed_file = "pipeline_test_data.bin.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate a binary file in non-fastq non-sam format
        binary_data = bytearray(random.getrandbits(8) for _ in range(1024 * 1024))  # 1MB random binary data
        self.original_hash = hashlib.md5(binary_data).hexdigest()
        
        with open(self.source_file, 'wb') as f:
            f.write(binary_data)
        
        
        # Add test command: pipe input, file output
        compress_command = f"cat {self.source_file} | ./release-release/bin/pbgz compress -o {self.compressed_file}"
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
        with open(self.source_file, 'rb') as f1, open(self.decompressed_file, 'rb') as f2:
            original_data = f1.read()
            decompressed_data = f2.read()

        # Use hash value comparison to avoid memory comparison issues
        decompressed_hash = hashlib.md5(decompressed_data).hexdigest()
        original_hash = hashlib.md5(original_data).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = BinaryPipelineTest()
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
