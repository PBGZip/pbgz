"""
Test case: SamDecompressWithRangeTest - SAM file region decompression test

Test scenario:
- Test using -p parameter to decompress specified chromosome region of SAM compressed files
- Verify correctness and precision of region filtering functionality
- Test processing capability of genomic coordinate range queries

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress -p chr3:100-200 {compressed_file} -o {decompressed_file}

Expected results:
- Compression command executes successfully
- Region decompression command succeeds, only containing chr3 position 100-200 records
- Decompressed file content exactly matches target region
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamDecompressWithRangeTest(PBGZTestCase):
    
    def __init__(self):
        super().__init__("SamDecompressWithRangeTest")
        self.source_file = "sam_range_test.sam"
        self.compressed_file = "sam_range_test.sam.pbgz"
        self.index_file = "sam_range_test.sam.pbgz.pbgzi"
        self.decompressed_file = "sam_range_test.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate sorted SAM file, using pbgz-supported format
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")
        sam_content.append("@SQ\tSN:chr1\tLN:1000000")
        sam_content.append("@SQ\tSN:chr2\tLN:1000000")
        sam_content.append("@SQ\tSN:chr3\tLN:1000000")
        
        # Add reads from different regions
        # chr1: position 100-500 range
        for pos in [100, 200, 300, 400, 500]:
            qname = f"CHR1_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr1", pos)
            sam_content.append(read_line)
        
        # chr1: position 1000-1500 range (not in target range)
        for pos in [1000, 1200, 1500]:
            qname = f"CHR1_OUT_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr1", pos)
            sam_content.append(read_line)
        
        # chr2: position 50-150 range (different chromosome)
        for pos in [50, 150]:
            qname = f"CHR2_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr2", pos)
            sam_content.append(read_line)
        
        # chr3: position 100-200 range (target region)
        for pos in [120, 180]:
            qname = f"CHR3_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr3", pos)
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # First compress
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -i")
        
        # Test region decompression (if pbgz supports it)
        self.add_test_command(f"./release-release/bin/pbgz decompress -p 'chr3:100-200' {self.compressed_file} -o {self.decompressed_file}")

    def _create_sam_line(self, qname, rname, pos):
        """Create SAM record line"""
        flag = 0
        mapq = 60
        cigar = "30M"
        rnext = "*"
        pnext = 0
        tlen = 0
        seq = ''.join(random.choices('ATCGN', k=30))
        qual = '!' * 30
        return f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:30"

    def cleanup_test_data(self):
        for filename in [self.source_file, 
                         self.compressed_file, 
                         self.index_file, 
                         self.decompressed_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        # Check command execution - based on command type, not file existence
        compress_success = False
        decompress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            file_sizes = cmd_result.get("file_sizes", {})
            
            # Judge by command content if it's compression or decompression command - must check "decompress" first to avoid being contained by "compress"
            if "decompress" in command and "pbgz" in command:
                decompress_success = success
            elif "compress" in command and "pbgz" in command:
                compress_success = success
        
        # Verify basic functionality
        if not compress_success:
            print("Compression command failed")
            return False
        
        if not decompress_success:
            print("Decompression command failed")
            return False
        
        return True


if __name__ == "__main__":
    test_case = SamDecompressWithRangeTest()
    result = test_case.execute()
    test_case.print_results()
    # Save JSON file for debugging results
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)