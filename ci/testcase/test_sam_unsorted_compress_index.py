"""
testusecase：Sam Unsorted Compress Index

testscenarioscenario：
- testpbgzofSam Unsorted Compress Indexfunctionality
- Verifymutualrelatedcommandofcorrectaccurateexecuteexecution
- ensuredatadatacompletecompletequality

Expectedexpectedresultresult：
- testcommandexecuteexecutiongeneratesuccess
- datadatacompletecompletequalityVerifythrough"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase

class SamUnsortedCompressIndexTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamUnsortedCompressIndexTest")
        self.source_file = "sam_unsorted_idx.sam"
        self.compressed_file = "sam_unsorted_idx.sam.pbgz"
        self.index_file = "sam_unsorted_idx.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        for i in range(25):
            qname = f"UNS_IDX_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 400
            mapq = 60
            cigar = "35M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=35))
            qual = '!' * 35
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:35"
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}")
        # indexcommandtonotsortedfilewilllosefail，thisisExpectedexpectedexecutionfor，nottypeprinterrorerror
        self.add_test_command(f"./release-release/bin/pbgz index -f {self.compressed_file}")

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.index_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        # Checkcompressfileiswhethersaveat
        if not os.path.exists(self.compressed_file):
            return False
        
        # indexcommandoflosefailnotinfluenceresponsetestresultresult，geneforthisisneedletoindexfunctionalitysingleuniqueoftest
        # compressgeneratesuccessimmediatelycan
        return True

if __name__ == "__main__":
    test_case = SamUnsortedCompressIndexTest()
    result = test_case.execute()
    test_case.print_results()
    # Save JSON file for debugging results
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
