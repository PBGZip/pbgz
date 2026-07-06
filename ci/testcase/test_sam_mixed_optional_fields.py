"""
testtestusecase：SamMixedOptionalFieldsTest - mixmatchcanselectfieldfieldSAMfilefilecompresscompressdecompresscompresstesttest

testtestscenarioscenario：
- testtestSAMfilefileincontainingcontainmixmatchcanselectfieldfield（headerdividerecordrecordhave，headerdividerecordrecordnohave）
- Verifyverifypbgztocanselectfieldfieldhandlemanageofonecausequality
- testtestdatadataformatformatflexibleactivequality

useusecommand：
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Expectedexpectedresultresult：
- compresscompressgeneratesuccess，correctaccuratehandlemanagemixmatchcanselectfieldfield
- decompresscompressgeneratesuccess，canselectfieldfieldcompletecompleteensurekeep
- MD5verifyVerifycompletefullonecause
"""

import sys
import random
import hashlib

# AddparenttargetrecordtoPythonroadpath
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase
import random
import hashlib
# AddparenttargetrecordtoPythonroadpath
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase
# haveheavycomplexsys.path.insert，needwantrepaircomplexonedownresultstructure，butfirstcontinuecontinuehandlemanageitsotherfilefile
# AddparenttargetrecordtoPythonroadpath
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class SamMixedOptionalFieldsTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamMixedOptionalFieldsTest")
        self.source_file = "sam_mixed_opt.sam"
        self.compressed_file = "sam_mixed_opt.sam.pbgz"
        self.decompressed_file = "sam_mixed_opt.sam.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        for i in range(20):
            qname = f"MIXED_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 200
            mapq = 60
            cigar = "25M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=25))
            qual = '!' * 25
            # mixmatchcanselectfieldfield
            if i % 3 == 0:
                read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            elif i % 3 == 1:
                read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0"
            else:
                read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:25\tAS:i:25"
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
    test_case = SamMixedOptionalFieldsTest()
    result = test_case.execute()
    test_case.print_results()
    # foradjusttestresultresultensuresaveJSONfilefile
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)
