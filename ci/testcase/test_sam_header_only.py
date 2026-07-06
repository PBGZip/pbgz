"""
Test case: SamHeaderOnlyTest - Header-only SAM file compression/decompression test

Test scenario:
- Test SAM file processing containing only header, no actual records
- Verify pbgz's capability to handle SAM files with empty data body
- Test file boundary cases and minimum dataset processing

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Parameter description:
- SAM header lines starting with @ for metadata
- Test extreme case with only header, no sequence records
- Verify compression/decompression stability of minimum dataset

Expected results:
- Compression succeeds, processing header-only file
- Decompression succeeds, header completely preserved
- MD5 verification completely matches
"""

import os
import sys
import hashlib
from testcase.pbgz_test_framework import PBGZTestCase

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class SamHeaderOnlyTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamHeaderOnlyTest")
        self.source_file = "sam_header_only.sam"
        self.compressed_file = "sam_header_only.sam.pbgz"
        self.decompressed_file = "sam_header_only.sam.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate header-only SAM file
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ref1\tLN:1000")
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}")
        self.add_test_command(f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}")
    
    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.decompressed_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print(f"Warning: Failed to remove {filename}: {e}")
    
    def verify_expected_results(self) -> bool:
        if not os.path.exists(self.compressed_file):
            print(f"Error: Compressed file {self.compressed_file} not created")
            return False
        
        if not os.path.exists(self.decompressed_file):
            print(f"Error: Decompressed file {self.decompressed_file} not created")
            return False
        
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        if original_content == decompressed_content:
            
            return True
        else:
            print(f"✗ Header-only SAM file test failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = SamHeaderOnlyTest()
    result = test_case.execute()
    test_case.print_results()
    # Save JSON file for debugging results
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)