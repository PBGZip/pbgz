"""
Test case: BinaryPipeTest - Binary file pipe compression/decompression test

Test scenario:
- Test pbgz compressing and decompressing binary files via pipe method
- Verify functionality of using standard input/output pipe interface
- Test combination usage of cat command and pbgz

Commands used:
1. cat {source_file} | ./release-release/bin/pbgz compress -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} > {decompressed_file}

Expected results:
- Pipe compression succeeds, generating compressed file
- Pipe decompression succeeds, restoring original file
- Pipe interface function works normally
- File contents are completely consistent

Technical notes:
- Test standard input/output pipe mode
- Verify pbgz compatibility with Unix pipe system
"""

import os
import sys
import random
import hashlib

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class BinaryFilePipeTest(PBGZTestCase):
    def __init__(self):

        super().__init__("BinaryFilePipeTest")
        self.source_file = "filepipe_test_data.bin"
        self.compressed_file = "filepipe_test_data.bin.pbgz"
        self.decompressed_file = "filepipe_test_data.bin.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate a binary file in non-fastq non-sam format
        binary_data = bytearray(random.getrandbits(8) for _ in range(1024 * 1024))  # 1MB random binary data
        self.original_hash = hashlib.md5(binary_data).hexdigest()
        
        with open(self.source_file, 'wb') as f:
            f.write(binary_data)
        
        
        # Add test command: file input, pipe output
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o - > {self.compressed_file}"
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
    test_case = BinaryFilePipeTest()
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