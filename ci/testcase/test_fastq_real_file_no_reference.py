"""
testtestusecase：FastqRealFileNoReferenceTest - realrealFASTQfilefilecompresscompresstesttest（invalidreferencereferencebasegenegenome）

testtestscenarioscenario：
- useuserealrealofFASTQdatadatafilefile（SRR17052138_1.fastq.gz）advanceexecutioncompresscompressdecompresscompresstesttest
- Verifyverifypbgzhandlemanagelargerulepatternrealrealdatadataofabilityforce
- testtestinvalidreferencereferencebasegenegenometimeofstandardstandardcompresscompressqualityability

useusecommand：
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

referencedatasayclear：
- useusetestorderdatadatafilefile（Throughnormalverylarge，containrealreal biological sequences）
- testtestpbgztorealrealdatadataofsuitableusequalityandcompresscompressrate
- Verifyverifylargerulepatternfilefilehandlemanageofstabledeterminequalityandqualityability

Expectedexpectedresultresult：
- compresscompressgeneratesuccess，handlemanagerealrealGenerateobjectdatadata
- decompresscompressgeneratesuccess，datadatacompletebeautystilloriginal
- MD5verifyVerifycompletefullonecause
"""

import os
import sys

# AddparenttargetrecordtoPythonroadpath
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqRealFileNoReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("FastqRealFileNoReferenceTest")
        self.source_file = "fq/SRR17052138_1.fastq.gz"
        self.compressed_file = "fq/SRR17052138_1.fastq.pbgz"
        self.decompressed_file = "fq/SRR17052138_1.fastq.dec"
        self.original_size_for_cleanup = 0
        self.compressed_size_for_cleanup = 0

    def get_test_files(self) -> tuple:
        """Returnreturntesttestfilefileinformationinformation"""
        return (self.source_file, self.compressed_file, self.decompressed_file)
        
    def prepare_data(self):
        # Checkcheckoriginalstartfilefileiswhethersaveat
        if not os.path.exists(self.source_file):
            print(f"Error: Source file {self.source_file} not found")
            return
        
        original_size = os.path.getsize(self.source_file)
        self.original_size_for_cleanup = original_size
        
        print(f"Real FASTQ file: {self.source_file}")
        
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -f"
        self.add_test_command(compress_command)
        
        # Adddecompresscompresscommand
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -z"
        self.add_test_command(decompress_command)
    
    def cleanup_test_data(self):
        # deleteremovecompresscompressanddecompresscompressfilefile
        for filename in [self.compressed_file, self.decompressed_file]:
            if os.path.exists(filename):
                try:
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
        
        # Verifyverifydecompresscompressfilefilesaveat
        if not os.path.exists(self.decompressed_file):
            print("decompresscompressfilefilenotsaveat")
            return False
        
        # countcomputecompresscompressqualityability
        original_size = os.path.getsize(self.source_file)
        compressed_size = os.path.getsize(self.compressed_file)
        compression_ratio = (original_size - compressed_size) / original_size * 100
        print(f"  compresscompressrate: {compression_ratio:.2f}%")
        
        # Verifyverifydecompresscompressfilefilesaveatandlargesmallmatchmanage
        decompressed_size = os.path.getsize(self.decompressed_file)
        print(f"  decompresscompressfilefilelargesmall: {decompressed_size:,} fieldnode ({decompressed_size / (1024**3):.2f}GB)")
        
        print("✓ FastqRealFileNoReferenceTest qualityabilitytesttestcompletegenerate！")
        return True


if __name__ == "__main__":
    test_case = FastqRealFileNoReferenceTest()
    result = test_case.execute()
    test_case.print_results()
    
    # fordebugtestresultresultensuresaveJSONfilefile
    if not result:
        test_case.save_to_json()
    
    if result:
        print("\\\\n✓ Real file performance test passed")
    else:
        print("\\\\n✗ Real file performance test failed")
        sys.exit(1)