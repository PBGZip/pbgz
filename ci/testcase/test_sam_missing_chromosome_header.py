"""
testtestusecase：SamMissingChromosomeInHeaderTest - SAMfilefilelackfewchromosomechromosomebodyheaderheaderdeterminemeaningtesttest

testtestscenarioscenario：
- testtestSAMfilefileinrealactualdatadatacontainingcontainofchromosomechromosomebodynohaveat@SQ headerindeterminemeaningofsituationsituation
- VerifyverifypbgztonotcompletecompleteSAMformatformatofcapacityerrorabilityforce
- testtestchromosomechromosomebodyindexusecompletecompletequalityCheckcheckofedgeboundaryhandlemanage

useusecommand：
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Expectedexpectedresultresult：
- compresscompressgeneratesuccess，pbgzabilityenoughhandlemanagelackfewchromosomechromosomebodydeterminemeaningofdatadata
- decompresscompressgeneratesuccess，placehavedatadatacompletecompleteensurekeep
- MD5verifyVerifycompletefullonecause
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamMissingChromosomeInHeaderTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamMissingChromosomeInHeaderTest")
        self.source_file = "sam_missing_chr.sam"
        self.compressed_file = "sam_missing_chr.sam.pbgz"
        self.decompressed_file = "sam_missing_chr.sam.dec"

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        sam_content = []
        
        # Headerinonlydeterminemeaningheaderdividechromosomechromosomebody
        sam_content.append("@HD\tVN:1.0\tSO:coordinate")
        sam_content.append("@SQ\tSN:chr1\tLN:1000000")  # determinemeaningpast tense markerchr1
        sam_content.append("@SQ\tSN:chr2\tLN:500000")   # determinemeaningpast tense markerchr2
        sam_content.append("@SQ\tSN:chr3\tLN:2000000")  # determinemeaningpast tense markerchr3
        # notemeaning：thisplacenohavedeterminemeaningchr4
        
        # Addcontainingcontainalreadydeterminemeaningchromosomechromosomebodyofdatadata
        for pos in [100, 200, 300]:
            qname = f"CHR1_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr1", pos)
            sam_content.append(read_line)
        
        for pos in [50, 150]:
            qname = f"CHR2_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr2", pos)
            sam_content.append(read_line)
        
        # relatedkeypoint：Addcontainingcontainnotatheaderindeterminemeaningofchromosomechromosomebodyofdatadata（chr4）
        for pos in [1000, 1200]:
            qname = f"CHR4_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr4", pos)  # chr4notatheaderindeterminemeaning
            sam_content.append(read_line)
        
        # againAddonesomealreadydeterminemeaningchromosomechromosomebodyofdatadata
        for pos in [400, 500]:
            qname = f"CHR3_{pos:04d}"
            read_line = self._create_sam_line(qname, "chr3", pos)
            sam_content.append(read_line)
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # compresscompress
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}")
        
        # decompresscompress
        self.add_test_command(f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}")

    def _create_sam_line(self, qname, rname, pos):
        """CreatecreateSAMrecordrecordexecution"""
        flag = 0
        mapq = 60
        cigar = "30M"
        rnext = "*"
        pnext = 0
        tlen = 0
        seq = ''.join(random.choices('ATCGN', k=30))
        qual = '!' * 30
        return f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:30"

    def cleanup_test_data(self):
        for filename in [self.source_file, self.compressed_file, self.decompressed_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception:
                pass

    def verify_expected_results(self) -> bool:
        # Checkcheckcommandexecuteexecutionsituationsituation
        compress_success = False
        decompress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            # baseatcommandclasstypepreciseaccuratejudgebreak
            parts = command.split()
            if len(parts) >= 2:
                main_cmd = parts[1]
                
                if main_cmd == "compress":
                    compress_success = success
                elif main_cmd == "decompress":
                    decompress_success = success
        
        # Verifyverifybasebasesuccessability
        if not compress_success or not decompress_success:
            print(f"compresscompressgeneratesuccess: {compress_success}, decompresscompressgeneratesuccess: {decompress_success}")
            return False
        
        # Verifyverifycompresscompressfilefilesaveat
        if not os.path.exists(self.compressed_file):
            print("compresscompressfilefilenotsaveat")
            return False
        
        # Verifyverifydecompresscompressfilefilesaveat
        if not os.path.exists(self.decompressed_file):
            print("decompresscompressfilefilenotsaveat")
            return False
        
        # Verifyverifydatadatacompletecompletequality
        with open(self.source_file, 'r') as f:
            original_lines = f.readlines()
        
        with open(self.decompressed_file, 'r') as f:
            decompressed_lines = f.readlines()
        
        if len(original_lines) != len(decompressed_lines):
            print(f"executiondatanotmatchpair: originalstart{len(original_lines)}execution vs decompresscompress{len(decompressed_lines)}execution")
            return False
        
        # MD5verifyVerify
        with open(self.source_file, 'rb') as f:
            original_md5 = hashlib.md5(f.read()).hexdigest()
        
        with open(self.decompressed_file, 'rb') as f:
            decompressed_md5 = hashlib.md5(f.read()).hexdigest()
        
        if original_md5 != decompressed_md5:
            print(f"MD5notmatchpair: {original_md5} vs {decompressed_md5}")
            return False
        
        # Verifyverifyrelatedkeypoint：containingcontainnotdeterminemeaningchromosomechromosomebodyofrecordrecord
        has_un_defined_chr = False
        for line in decompressed_lines:
            if line.startswith('@'):  # skipthroughheader
                continue
            fields = line.strip().split('\t')
            if len(fields) >= 3 and fields[2] == "chr4":
                has_un_defined_chr = True
                break
        
        if not has_un_defined_chr:
            print("decompresscompressfilefileinlackfewnotdeterminemeaningchromosomechromosomebodyofrecordrecord(chr4)")
            return False
        
        print("SAMfilefilecontainingcontainnotatheaderindeterminemeaningofchromosomechromosomebodytimealsoabilitycorrectnormalcompresscompressdecompresscompress")
        return True


if __name__ == "__main__":
    test_case = SamMissingChromosomeInHeaderTest()
    result = test_case.execute()
    test_case.print_results()
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)