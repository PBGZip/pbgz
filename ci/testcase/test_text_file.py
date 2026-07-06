"""
Test case: TextFileTest - Text file compression/decompression test

Test scenario:
- Test pbgz's compression and decompression for normal text files
- Verify integrity preservation of text content
- Test handling of diverse characters (letters, numbers, symbols)

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Parameter description:
- Test file contains 1000 lines of text, each line contains diverse characters
- Characters include: uppercase letters, lowercase letters, numbers, special symbols
- Each line length is approximately 100 characters

Expected results:
- Compression succeeds, generating compressed file
- Decompression succeeds, restoring original text file
- Decompressed file MD5 matches original file completely
- All characters and line structure correctly preserved
- Compression ratio is reasonable (because text has repeated patterns)

Technical notes:
- testmostnormaluseoffilebasefileclasstype
- VerifypbgzhandlingThroughusefilebasedatadataofabilityforce
- filebasefilecompressisbasebasefunctionalityVerify
"""

import os
import sys
import random
import hashlib

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class TextFileTest(PBGZTestCase):
    def __init__(self):

        super().__init__("TextFileTest")
        self.source_file = "text_file_data.txt"
        self.compressed_file = "text_file_data.txt.pbgz"
        self.decompressed_file = "text_file_data.txt.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generatefilebasefileinsidecapacity
        full_text = []
        num_lines = 1000  # Generate 1000 lines of text
        
        for i in range(num_lines):
            line = f"This is line {i} of text file data for compression testing. "
            line += "Each line contains various characters: ABCDEF123456!@#$%^&*()"
            full_text.append(line)
        
        self.original_hash = hashlib.md5('\n'.join(full_text).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in full_text:
                f.write(line + '\n')
        
        
        # Addcompresstestcommand
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command, cmd_id=1)
        
        # Adddecompresstestcommand
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command, cmd_id=2)
    
    def cleanup_test_data(self):
        """clearmanagetestusecaseCreateoftemporarytimefile"""
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
        
        decompressed_hash = hashlib.md5(decompressed_content).hexdigest()
        original_hash = hashlib.md5(original_content).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = TextFileTest()
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