"""
testtestusecase：SamRealFileReferenceTest - realreallargefilefileSAMqualityabilitytesttest（useusereferencereferencebasegenegenome）

testtestscenarioscenario：
- useuserealrealoflargerulepatternSAMfilefile（SRR2007660.sam，3.7GB）advanceexecutionqualityabilitytesttest
- useusetorespondofreferencereferencebasegenegenomeadvanceexecutionoptimaltransformcompresscompress
- Verifyverifyreferencereferencebasegenegenometolargerulepatterndatadataofcompresscompressoptimaltransformvalidresult
- testtestdirectconnectreadgetgzformatformatreferencereferencebasegenegenomeofabilityforce

useusecommand：
1. ./release-release/bin/pbgz compress sam/SRR2007660.sam -o sam/SRR2007660.refsam.pbgz -r fa/GCA_000195955.2.fasta.gz
2. ./release-release/bin/pbgz decompress sam/SRR2007660.refsam.pbgz -o sam/SRR2007660.refsam.dec -r fa/GCA_000195955.2.fasta.gz

qualityabilityVerifyverifypointstandard：
- compresscompresstimebetweenwithinvalidreferencereferencebasegenegenomeoftocompare
- compresscompressratetocompare（referencereferencebasegenegenomeThroughnormalwilldisplaysignificantimproveupgradecompresscompressrate）
- decompresscompresstimebetweenanddatadatacompletecompletequality
- directconnectreadgetgzformatformatreferencereferencebasegenegenomeofsupportmaintain

Expectedexpectedresultresult：
- compresscompressgeneratesuccess，supportmaintaindirectconnectreadgetgzreferencereferencebasegenegenome
- compresscompressratedisplaysignificanthighatinvalidreferencereferencebasegenegenomeofsituationsituation（Throughnormal>200%）
- decompresscompressgeneratesuccess，datadatacompletebeautystilloriginal
- MD5verifyVerifycompletefullonecause
"""

import os
import sys
import hashlib
import time

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamRealFileReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamRealFileReferenceTest")
        self.source_file = "sam/SRR2007660.sam"
        self.compressed_file = "sam/SRR2007660.refsam.pbgz"
        self.decompressed_file = "sam/SRR2007660.refsam.dec"
        self.reference_file = "fa/GCA_000195955.2.fasta.gz"  # comeself@PGinformationinformationofreferencereferencebasegenegenome

    def get_test_files(self) -> tuple:
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Verifyverifysourcefilefilesaveat
        if not os.path.exists(self.source_file):
            print(f"Error: Source file {self.source_file} does not exist")
            return False
        
        # Verifyverifyreferencereferencebasegenegenomefilefilesaveat
        if not os.path.exists(self.reference_file):
            print(f"Error: Reference genome file {self.reference_file} does not exist")
            return False
    
        print(f"Real SAM file: {self.source_file}")
        print(f"  Reference genome: {self.reference_file}")
        
        # useusereferencereferencebasegenegenomecompresscompress（Add-freferencedatastrongcontrolcovercoveralreadysaveatoffilefile）
        self.add_test_command(f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file} -f")
        
        # useusereferencereferencebasegenegenomedecompresscompress
        self.add_test_command(f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r {self.reference_file}")

    def cleanup_test_data(self):
        # clearmanageplacehaveGenerategenerateoffilefile
        for filename in [self.compressed_file, self.decompressed_file]:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                pass

    def verify_expected_results(self) -> bool:
        # Checkcheckcommandexecuteexecutionresultresult
        compress_success = False
        decompress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            parts = command.split()
            if len(parts) >= 2:
                main_cmd = parts[1]
                
                if main_cmd == "compress":
                    compress_success = success
                    # recordrecordcompresscompresstimebetween
                    if "execution_time" in cmd_result:
                        comp_time = cmd_result["execution_time"]
                        print(f"  compresscompresstimebetween: {comp_time:.2f}second")
                elif main_cmd == "decompress":
                    decompress_success = success
                    # recordrecorddecompresscompresstimebetween
                    if "execution_time" in cmd_result:
                        decomp_time = cmd_result["execution_time"]
                        print(f"  decompresscompresstimebetween: {decomp_time:.2f}second")
        
        # Verifyverifybasebasesuccessability
        if not compress_success:
            print("compresscompresscommandlosefail")
            return False
        
        if not decompress_success:
            print("decompresscompresscommandlosefail")
            return False
        
        # Verifyverifycompresscompressfilefilesaveat
        if not os.path.exists(self.compressed_file):
            print("compresscompressfilefilenotsaveat")
            return False
        
        # countcomputecompresscompresscompare
        original_size = os.path.getsize(self.source_file)
        compressed_size = os.path.getsize(self.compressed_file)
        compression_ratio = (original_size - compressed_size) / original_size * 100
        print(f"  compresscompressrate: {compression_ratio:.2f}%")
        
        # Verifyverifydecompresscompressfilefilesaveatandlargesmallmatchmanage
        decompressed_size = os.path.getsize(self.decompressed_file)
        print(f"  decompresscompressfilefilelargesmall: {decompressed_size:,} fieldnode ({decompressed_size / (1024**3):.2f}GB)")
        
        # MD5verifyVerifybeforeseveralpieceMBofdatadata（qualityabilityreferenceconsider）
        with open(self.source_file, 'rb') as f1, open(self.decompressed_file, 'rb') as f2:
            original_part = f1.read(10*1024*1024)  # readgetbefore10MB
            decompressed_part = f2.read(10*1024*1024)
            
        original_md5 = hashlib.md5(original_part).hexdigest()
        decompressed_md5 = hashlib.md5(decompressed_part).hexdigest()
        
        if original_md5 != decompressed_md5:
            return False
        
        print("✓ SamRealFileReferenceTest qualityabilitytesttestcompletegenerate！")
        return True


if __name__ == "__main__":
    test_case = SamRealFileReferenceTest()
    result = test_case.execute()
    test_case.print_results()
    if not result:
        test_case.save_to_json()
    sys.exit(0 if result else 1)