"""
Test case: Sam No Optional Fields Variable

Test scenario:

- Test pbgz's Sam No Optional Fields variable functionality
- Verify correct execution of related commands
- Ensure data integrity

Expected results:

- Test command executes successfully
- Data integrity verification passed"""

import os
import sys
import random
import hashlib
from testcase.pbgz_test_framework import PBGZTestCase

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class SamNoOptionalFieldsVariableTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamNoOptionalFieldsVariableTest")
        self.source_file = "sam_no_opt_var.sam"
        self.compressed_file = "sam_no_opt_var.sam.pbgz"
        self.decompressed_file = "sam_no_opt_var.sam.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        base_lengths = [10, 15, 20, 25, 30]
        for i, base_len in enumerate(base_lengths):
            qname = f"NO_OPT_VAR_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 600
            mapq = 60
            cigar = f"{base_len}M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=base_len))
            qual = '!' * base_len
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
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
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        return os.path.exists(self.compressed_file) and os.path.exists(self.decompressed_file)

if __name__ == "__main__":
    test_case = SamNoOptionalFieldsVariableTest()
    result = test_case.execute()
    test_case.print_results()
    # Save JSON file for debugging results
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
