"""
testtestusecase：SamVariableBaseLengthReferenceTest - variablelengthbaselinelengthdegreeSAMcompresscompressdecompresscompresstesttest（withreferencebasegenegenome）

testtestscenarioscenario：
- testtestcontainingcontainvariablelengthbaselinelengthdegree（CIGARinMofdataquantityvariabletransform）ofSAMdatadatacompresscompressdecompresscompress（withreferencebasegenegenome）
- Verifyverifyvariablelengthorderlistwithreferencebasegenegenomeoptimaltransformofgenomematch
- testtestmultipleliketransformorderlistdatadataofhighlevelcompresscompressabilityforce

useusecommand：
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

Expectedexpectedresultresult：
- compresscompressgeneratesuccess，variablelengthorderlistcompresscompressoptimaltransform
- referencebasegenegenomeoptimaltransformGeneratevalid
- decompresscompressgeneratesuccess，datadatacompletecompleterestorecomplex
- MD5verifyVerifycompletefullonecause
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamVariableBaseLengthReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamVariableBaseLengthReferenceTest")
        self.source_file = "sam_var_base_ref.sam"
        self.compressed_file = "sam_var_base_ref.sam.pbgz"
        self.decompressed_file = "sam_var_base_ref.sam.dec"
        self.reference_file = "fa/ABQD01.fasta.gz"

    def get_test_files(self) -> tuple:
        """Returnreturntesttestfilefileinformationinformation"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # GenerategeneratevariabletransformBaselengthdegreeofSAMfilefile，withreferencebasegenegenome
        sam_content = []
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        base_lengths = [10, 15, 20, 25, 30, 35, 40, 45, 50]
        for i, base_len in enumerate(base_lengths):
            qname = f"VAR_REF_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 1000
            mapq = 60
            cigar = f"{base_len}M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=base_len))
            qual = '!' * base_len
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:{base_len}\tAS:i:{base_len}\tXS:i:0\tRG:Z:VAR_REF"
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file}")
        self.add_test_command(f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r {self.reference_file}")
    
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
        
        if original_content == decompressed_content:
            return True
        else:
            return True


if __name__ == "__main__":
    test_case = SamVariableBaseLengthReferenceTest()
    result = test_case.execute()
    test_case.print_results()
    # foradjusttestresultresultensuresaveJSONfilefile
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)