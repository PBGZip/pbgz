"""
Test case: SamCompressWithIndexVariableTest - Index generation test for ordered variable-length SAM files

Test scenario:
- Test that ordered SAM files containing variable-length sequences can successfully generate index files
- Verify coordination between ordered data and index creation
- Test index processing for variable-length sequences

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz index -f {compressed_file}

Expected results:
- Compression succeeds
- Index file successfully generated and not empty
- Index file format is correct (5 tab-separated fields per line)
"""

import sys
import random
import hashlib

import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase

class SamCompressWithIndexVariableTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamCompressWithIndexVariableTest")
        self.source_file = "sam_index_var.sam"
        self.compressed_file = "sam_index_var.sam.pbgz"
        self.index_file = "sam_index_var.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")  # Change to ordered
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        base_lengths = [15, 20, 25, 30, 35]
        for i, base_len in enumerate(base_lengths):
            qname = f"INDEX_VAR_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 700  # Position increasing, ensuring ordered
            mapq = 60
            cigar = f"{base_len}M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=base_len))
            qual = '!' * base_len
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:{base_len}"
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}")
        # Ordered files should be able to successfully generate index files
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
            
            # Precise judgment based on command type, avoiding string inclusion issues
            parts = command.split()
            if len(parts) >= 2:
                main_cmd = parts[1]
                
                if main_cmd == "index":
                    index_success = success
                elif main_cmd == "compress":
                    compress_success = success
        
        # Compression must succeed
        if not compress_success:
            print("Compression command failed")
            return False
        
        # Compressed file must exist
        if not os.path.exists(self.compressed_file):
            print(f"Compressed file {self.compressed_file} not found")
            return False
        
        # For ordered SAM files, we expect index command to succeed
        if not index_success:
            print("Index command failed - ordered SAM should be able to generate index")
            return False
        
        # Index file must exist and not be empty
        if not os.path.exists(self.index_file):
            print(f"Index file {self.index_file} not created")
            return False
        
        if os.path.getsize(self.index_file) == 0:
            print(f"Index file {self.index_file} is empty")
            return False
        
        # Verify index file format: 5 tab-separated fields per line
        try:
            with open(self.index_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    fields = line.split('\t')
                    if len(fields) != 5:
                        print(f"Index file line {line_num} has incorrect field count: expected 5, got {len(fields)}")
                        print(f"  Content: {line}")
                        return False
        except Exception as e:
            print(f"Error reading index file: {e}")
            return False
        
        return True

if __name__ == "__main__":
    test_case = SamCompressWithIndexVariableTest()
    result = test_case.execute()
    test_case.print_results()
    # Save JSON file for debugging results
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
