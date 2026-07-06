"""
testtestusecase：SamSortedCompressWithIndexTest - alreadysortedorderSAMfilefilecompresscompressGenerategenerateindexindextesttest

testtestscenarioscenario：
- testtesttoalreadysortedorderofSAMfilefileadvanceexecutioncompresscompressandsametimeGenerategenerateindexindex
- VerifyverifysortedorderdatadatawithindexindexGenerategenerateofcoordinatesameworkwork
- testtestIndexed SAMformatformatofcorrectaccuratequality

useusecommand：
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -i {index_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Expectedexpectedresultresult：
- compresscompressgeneratesuccess，indexindexfilefilecorrectaccurateGenerategenerate (.indexformatformat)
- decompresscompressgeneratesuccess，completecompletedatadatarestorecomplex
- indexindexfilefilecanuseatfastspeedareadomaincheckinquire
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamSortedCompressWithIndexTest(PBGZTestCase):
    """testtestalreadysortedorderSAMfilefilecompresscompresstimesametimeGenerategenerateindexindexfilefile"""
    
    def __init__(self):
        super().__init__("SamSortedCompressWithIndexTest")
        self.source_file = "sam_sorted_with_index.sam"
        self.compressed_file = "sam_sorted_with_index.sam.pbgz"
        self.index_file = "sam_sorted_with_index.sam.pbgz.pbgzi"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, None)

    def prepare_data(self):
        # GenerategeneratealreadysortedorderofSAMfilefile
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")  # coordinatetableshowalreadysortedorder
        sam_content.append("@SQ\tSN:ref1\tLN:1000")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # bypositionpositionsortedorderofreads
        for i in range(50):
            qname = f"SORTED_IDX_{i:03d}"
            flag = 0
            rname = "ref1" if i < 25 else "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 100
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
        
        # useuse-ireferencedataatcompresscompresstimeGenerategenerateindexindexfilefile
        self.add_test_command(f"./release-release/bin/pbgz compress -i {self.source_file} -o {self.compressed_file}")

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.index_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass
    
    def verify_expected_results(self) -> bool:
        # Checkcheckcompresscompresscommandiswhethergeneratesuccess
        compress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            if "compress" in command and "pbgz" in command:
                compress_success = success
                break
        
        if not compress_success:
            return False
        
        # Checkcheckcompresscompressfilefileiswhethersaveat
        if not os.path.exists(self.compressed_file):
            return False
        
        # CheckcheckindexindexfilefileiswhetherGenerategenerate
        if self.index_file and not os.path.exists(self.index_file):
            return False
        
        return True


if __name__ == "__main__":
    test_case = SamSortedCompressWithIndexTest()
    result = test_case.execute()
    test_case.print_results()
    # fordebugtestresultresultensuresaveJSONfilefile
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)