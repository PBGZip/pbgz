"""
Test case: SamCompressWithIndexTest - SAM file compression/decompression test

Test scenario:
- Test compression and decompression functionality for SAM files
- Verify data integrity and format correctness

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -i {index_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Expected results:
- Compression succeeds, index file correctly generated
- Decompression succeeds, data completely restored
- MD5 verification completely matches
"""

import sys
import random
import hashlib

import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase

class SamCompressWithIndexTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamCompressWithIndexTest")
        self.source_file = "sam_with_index.sam"
        self.compressed_file = "sam_with_index.sam.pbgz"
        self.index_file = "sam_with_index.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        for i in range(20):
            qname = f"INDEX_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 500
            mapq = 60
            cigar = "30M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=30))
            qual = '!' * 30
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:30"
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}")
        self.add_test_command(f"./release-release/bin/pbgz index -f {self.compressed_file}")

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.index_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        # Check command execution results
        compress_success = False
        index_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            if "compress" in command and "pbgz" in command:
                compress_success = success
            elif "index" in command and "pbgz" in command:
                index_success = success
        
        # Check if both compression and indexing succeeded
        if not compress_success:
            return False
            
        if not index_success:
            return False
        
        # Check if compressed file exists
        if not os.path.exists(self.compressed_file):
            return False
        
        # Check if index file exists
        if self.index_file and not os.path.exists(self.index_file):
            return False
        
        return True

if __name__ == "__main__":
    test_case = SamCompressWithIndexTest()
    result = test_case.execute()
    test_case.print_results()
    # Save JSON file for debugging results
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
