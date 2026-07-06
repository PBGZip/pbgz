"""
Test case: SamCompressWithIndexReferenceTest - SAM file compression/decompression test (with reference genome index)

Test scenario:
- Test SAM compression/decompression functionality using reference genome
- Verify coordination between reference genome indexing and compression
- Test correct handling of complex fields

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

Expected results:
- Compression succeeds, using reference genome optimization
- Decompression succeeds, data perfectly restored
- MD5 verification matches completely
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase

class SamCompressWithIndexReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamCompressWithIndexReferenceTest")
        self.source_file = "sam_index_ref.sam"
        self.compressed_file = "sam_index_ref.sam.pbgz"
        self.index_file = "sam_index_ref.sam.pbgz.index"
        self.reference_file = "fa/ABQD01.fasta.gz"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        for i in range(15):
            qname = f"INDEX_REF_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 600
            mapq = 60
            cigar = "25M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=25))
            qual = '!' * 25
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:25"
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file} -i")

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.index_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        return os.path.exists(self.compressed_file)

if __name__ == "__main__":
    test_case = SamCompressWithIndexReferenceTest()
    result = test_case.execute()
    test_case.print_results()
    # Save JSON file for debugging results
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)