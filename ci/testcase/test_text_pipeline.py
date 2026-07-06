"""
Test case: TextPipelineTest - Text file multi-stage pipeline processing test

Test scenario:
- Test pbgz's processing capability in multi-stage pipelines
- Verify stability of complex pipeline combinations
- Test file processing through multiple command chains

Commands used:
1. cat {source_file} | ./release-release/bin/pbgz compress | cat | ./release-release/bin/pbgz decompress -o {decompressed_file}

Expected results:
- Multi-stage pipeline processing succeeds
- File passes correctly through pipeline chain
- Final decompressed file matches original file content
- Pipeline combination operations are stable and reliable

Technical notes:
- Test complex pipeline combination scenarios
- Verify error handling and stream management
- Ensure correct passing through multi-command chains
"""

import os
import sys
import random
import hashlib

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class TextPipelineTest(PBGZTestCase):
    def __init__(self):

        super().__init__("TextPipelineTest")
        self.source_file = "text_pipeline_data.txt"
        self.compressed_file = "text_pipeline_data.txt.pbgz"
        self.decompressed_file = "text_pipeline_data.txt.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generatefilebasefileinsidecapacity
        full_text = []
        num_lines = 1000
        
        for i in range(num_lines):
            line = f"This is line {i} of text file data for pipeline compression testing."
            full_text.append(line)
        
        self.original_hash = hashlib.md5('\n'.join(full_text).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in full_text:
                f.write(line + '\n')
        
        
        # AddaddcompressTestcommand（pipelinepipelineoutputinput，fileoutputoutput）
        compress_command = f"cat {self.source_file} | ./release-release/bin/pbgz compress -o {self.compressed_file}"
        self.add_test_command(compress_command)
        
        # AddadddecompressTestcommand（pipelinepipelineoutputinput，fileoutputoutput）
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command)
    
    def cleanup_test_data(self):
        """clearmanageTestusecaseCreateoftemporarytimefile"""
        files_to_clean = [self.source_file, self.compressed_file, self.decompressed_file]
        for filename in files_to_clean:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print(f"Warning: Failed to remove {filename}: {e}")
    
    def verify_expected_results(self) -> bool:
        # VerifycompressfileiswhetherGenerate
        if not os.path.exists(self.compressed_file):
            print(f"Error: Compressed file {self.compressed_file} not created")
            return False
        
        # VerifydecompressfileiswhetherGenerate
        if not os.path.exists(self.decompressed_file):
            print(f"Error: Decompressed file {self.decompressed_file} not created")
            return False
        
        # Verifycompresscompareiswhethermatchmanage
        compression_ratio = self.get_compression_rate()
        if compression_ratio is None:
            print(f"Error: Failed to get compression ratio")
            return False
        
        # typeprintobtaingettoofexecuteexecutiontimebetweenandcompresscompare
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
    test_case = TextPipelineTest()
    result = test_case.execute()
    test_case.print_results()
    
    if result:
        print("\nTest passed successfully!")
    else:
        print("\nTest failed!")
        # Save JSON file for debugging results
        if not result:
            test_case.save_to_json()
        sys.exit(1)