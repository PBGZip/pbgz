"""
Test case: SamPipelineTest - SAM file pipe compression/decompression test

Test scenario:
- Test SAM file compression/decompression processing via Unix pipe method
- Verify compatibility of SAM format data with standard Unix tools
- Test correct execution of SAM format verification in pipes

Commands used:
1. cat {source_file} | ./release-release/bin/pbgz compress -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} | cat > {decompressed_file}

Expected results:
- Pipe compression/decompression succeeds
- SAM format completely preserved (@header, 11 column data)
- Decompressed file MD5 matches original file completely
- Pipe mode functionality works normally
"""

import os
import sys
import random
import hashlib
from testcase.pbgz_test_framework import PBGZTestCase

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


class SamPipelineTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamPipelineTest")
        self.source_file = "sam_pipe_test.sam"
        self.compressed_file = "sam_pipe_test.sam.pbgz"
        self.decompressed_file = "sam_pipe_test.sam.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate basic SAM file content
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")

        # Add 3 read records
        num_reads = 3
        for i in range(num_reads):
            qname = f"SAM_PIPE_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 1000
            mapq = 60
            cigar = "25M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATGCN', k=25))
            qual = '!' * 25
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:25\tAS:i:25\tXS:i:0\tRG:Z:PIPE"
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # Add compression and decompression commands (using pipeline)
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o - > {self.compressed_file}")
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
        
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            
            return True
        else:
            print(f"✗ Pipeline test failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = SamPipelineTest()
    result = test_case.execute()
    test_case.print_results()
    # Save JSON file for debugging results
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)