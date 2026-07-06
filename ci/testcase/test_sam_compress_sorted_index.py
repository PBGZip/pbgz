"""
Test case: SamCompressSortedIndexTest - Index creation test for already sorted SAM files

Test scenario:
- Test creating index files for already sorted SAM files
- Verify correct generation and format of index files
- Test index processing for sorted files (SO:coordinate)

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz index -f {compressed_file}

Expected results:
- Compression command executes successfully
- Index command executes successfully
- Index file exists and is not empty
- Index file format is correct: each line tab-separated, 5 fields

Technical notes:
- Sorted SAM files can create indexes
- Index file format is tab-separated 5-field format
- -f parameter forces overwrite of existing index files
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamCompressSortedIndexTest(PBGZTestCase):
    """Test index creation for sorted SAM files"""
    
    def __init__(self):
        super().__init__("SamCompressSortedIndexTest")
        self.source_file = "sam_sorted.sam"
        self.compressed_file = "sam_sorted.sam.pbgz"
        self.index_file = "sam_sorted.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        # Generate sorted SAM file
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")  # Note this is coordinate (sorted)
        sam_content.append("@SQ\tSN:ref1\tLN:1000")
        
        # Reads sorted by position
        for i in range(10):
            qname = f"SORTED_IDX_{i:03d}"
            flag = 0
            rname = "ref1"
            pos = (i + 1) * 100  # Position increasing, sorted
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
        
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}", cmd_id=1)
        self.add_test_command(f"./release-release/bin/pbgz index -f {self.compressed_file}", cmd_id=2)

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.index_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        # Check command execution results
        compress_success = index_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            # Check "index" first to avoid being contained by "compress"
            if "index" in command and "pbgz" in command:
                index_success = success
            elif "compress" in command and "pbgz" in command:
                compress_success = success
        
        # Sorted SAM file: both compression and indexing should succeed
        if not compress_success:
            return False
            
        # Index should succeed (because it's sorted file)
        if not index_success:
            return False
        
        # Verify index file content is not empty
        if not os.path.exists(self.index_file):
            print(f"Index file {self.index_file} not generated")
            return False
        
        if os.path.getsize(self.index_file) == 0:
            print(f"index file {self.index_file} is empty")
            return False
        
        # Verify index file content format: each line tab-separated, 5 fields
        with open(self.index_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:  # skip empty lines
                    continue
                
                fields = line.split('\t')
                if len(fields) != 5:
                    print(f"Index file line {line_num} has incorrect field count: expected 5, actual {len(fields)}")
                    print(f"  Content: {line}")
                    return False
        
        return True

if __name__ == "__main__":
    test_case = SamCompressSortedIndexTest()
    result = test_case.execute()
    test_case.print_results()
    # Save JSON file for debugging results
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)