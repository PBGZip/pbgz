"""
Test case: SamSortVariableBaseLengthTest - variable-lengthbaselinelengthofSAMsortedTest

Test scenario：
- TestSAMdatadataincontainingvariable-lengthbaselinelength（CIGARfieldcharacterstringinMofdataquantityvariabletransform）ofsorted
- Verifypbgzhandlingvariable-lengthsequenceofsortedalityforce
- Testmultiplesampletransformsequencedatadataofhandling

Usecommand：
1. ./release-release/bin/pbgz sort {input_file} -o {sorted_file}
2. ./release-release/bin/pbgz compress {sorted_file} -o {compressed_file}

Expected results：
- sortedcommandgeneratefunctionexecuteexecution
- variable-lengthsequencesortedcorrectEnsure
- compressdecompressgeneratefunction，datadatacompletecomplete
"""

import os
import sys
import random
import hashlib
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase

class SamSortVariableBaseLengthTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamSortVariableBaseLengthTest")
        self.source_file = "sam_sort_var.sam"
        self.sorted_file = "sam_sort_var.sort.sam"

    def get_test_files(self) -> tuple:
        return (self.source_file, None, None)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        base_lengths = [10, 15, 20, 25, 30]
        for i, base_len in enumerate(base_lengths):
            qname = f"SORT_VAR_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 1500
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
        
        self.add_test_command(f"./release-release/bin/pbgz sort {self.source_file} -o {self.sorted_file} -f")

    def cleanup_test_data(self):
        for filename in [self.source_file, self.sorted_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        return os.path.exists(self.sorted_file)

if __name__ == "__main__":
    test_case = SamSortVariableBaseLengthTest()
    result = test_case.execute()
    test_case.print_results()
    # fordebugtestresultresultensuresaveJSONfile
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)