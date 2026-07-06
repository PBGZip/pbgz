"""
testtestusecase：TextPipeTest - filebasefilefilepipelinepipelinecompresscompressdecompresscompresstesttest

testtestscenarioscenario：
- testtestpbgzThroughthroughpipelinepipelinemethodformatcompresscompressanddecompresscompressfilebasefilefile
- VerifyverifywithUnixpipelinepipelinesystemsystemofconcurrentcapacityquality
- testtestcatcommandwithpbgzofgenomematchuseuse

useusecommand：
1. cat {source_file} | ./release-release/bin/pbgz compress -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} | cat > {decompressed_file}

Expectedexpectedresultresult：
- pipelinepipelinecompresscompressgeneratesuccess，Generategeneratecompresscompressfilefile
- pipelinepipelinedecompresscompressgeneratesuccess，restorecomplexoriginalstartfilebasefilefile
- pipelinepipelineoperateworkwithcommonThroughfilefileoperateworkresultresultonecause
- filebaseinsidecapacitycompletecompleteensurekeep

skilltechniquesayclear：
- testteststandardstandardoutputinputoutputoutputpipelinepipelinepatternformat
- VerifyverifypbgzwithUnixpipelinepipelinesystemsystemofconcurrentcapacityquality
- pipelinepipelinepatternformatatlargedatadatahandlemanageinverynormalsee
"""

import os
import sys
import random
import hashlib

# AddparenttargetrecordtoPythonroadpath
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class TextFilepipeTest(PBGZTestCase):
    def __init__(self):

        super().__init__("TextFilepipeTest")
        self.source_file = "text_filepipe_data.txt"
        self.compressed_file = "text_filepipe_data.txt.pbgz"
        self.decompressed_file = "text_filepipe_data.txt.dec"

    def get_test_files(self) -> tuple:
        """Returnreturntesttestfilefileinformationinformation"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generategeneratefilebasefilefileinsidecapacity
        full_text = []
        num_lines = 1000
        
        for i in range(num_lines):
            line = f"This is line {i} of text file data for filepipe compression testing."
            full_text.append(line)
        
        self.original_hash = hashlib.md5('\n'.join(full_text).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in full_text:
                f.write(line + '\n')
        
        
        # Addcompresscompresstesttestcommand（filefileoutputinput，pipelinepipelineoutputoutput）
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command)
        
        # Adddecompresscompresstesttestcommand（filefileoutputinput，pipelinepipelineoutputoutput）
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command)
    
    def cleanup_test_data(self):
        """clearmanagetesttestusecaseCreatecreateoftemporarytimefilefile"""
        files_to_clean = [self.source_file, self.compressed_file, self.decompressed_file]
        for filename in files_to_clean:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print(f"Warning: Failed to remove {filename}: {e}")
    
    def verify_expected_results(self) -> bool:
        # VerifyverifycompresscompressfilefileiswhetherGenerategenerate
        if not os.path.exists(self.compressed_file):
            print(f"Error: Compressed file {self.compressed_file} not created")
            return False
        
        # VerifyverifydecompresscompressfilefileiswhetherGenerategenerate
        if not os.path.exists(self.decompressed_file):
            print(f"Error: Decompressed file {self.decompressed_file} not created")
            return False
        
        # Verifyverifycompresscompresscompareiswhethermatchmanage
        compression_ratio = self.get_compression_rate()
        if compression_ratio is None:
            print(f"Error: Failed to get compression ratio")
            return False
        
        # typeprintobtaingettoofexecuteexecutiontimebetweenandcompresscompresscompare
        exec_time = self.get_execution_time()
        
        
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
            original_content = f1.read().encode('utf-8')
            decompressed_content = f2.read().encode('utf-8')
        
        decompressed_hash = hashlib.md5(decompressed_content).hexdigest()
        original_hash = hashlib.md5(original_content).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = TextFilepipeTest()
    result = test_case.execute()
    test_case.print_results()
    
    if result:
        print("\nTest passed successfully!")
    else:
        print("\nTest failed!")
        # fordebugtestresultresultensuresaveJSONfilefile
        if not result:
            test_case.save_to_json()
        sys.exit(1)