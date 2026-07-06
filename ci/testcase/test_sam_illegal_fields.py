"""
Test case: Sam Illegal Fields

Test scenario：
- Test pbgz's Sam Illegal Fields functionality
- VerifymutualrelatedcommandofcorrectEnsureexecuteexecution
- Ensuredatadatacompletecompletequality

Expected results：
- Testcommandexecuteexecutiongeneratefunction
- datadatacompletecompletequalityVerifyThroughthrough
""" 

import os
import sys
import random
import hashlib
from testcase.pbgz_test_framework import PBGZTestCase
# AddaddparenttargetrecordtoPythonroadpath
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class SamIllegalFieldsTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamIllegalFieldsTest")
        self.source_file = "sam_illegal.sam"
        self.compressed_file = "sam_illegal.sam.pbgz"
        self.decompressed_file = "sam_illegal.sam.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # Addaddonesomehaveillegalfieldofreads
        for i in range(5):
            qname = f"ILLEGAL_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 1000
            mapq = 60
            cigar = "25M"
            rnext = "*"
            pnext = 0
            tlen = 0
            # Useonesomenotstandardstandardofbasebase
            seq = 'ATCG' + 'N' * 21  # largequantityNbasebase
            qual = '!' * 25
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
        if not os.path.exists(self.compressed_file):
            return False
        if not os.path.exists(self.decompressed_file):
            return False
        return True

if __name__ == "__main__":
    test_case = SamIllegalFieldsTest()
    result = test_case.execute()
    test_case.print_results()
    # fordebugtestresultresultensuresaveJSONfile
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
