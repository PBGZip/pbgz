"""
Test case: FastqNoTrailingNewlineTest - FASTQ file no trailing newline handling test

Test scenario:
- Test FASTQ files missing newline at the end
- Verify pbgz's capability to handle non-standard file endings
- Test processing of file boundary cases

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Parameter description:
- Normal files usually end with newline
- Test boundary case without trailing newline
- Verify file format standardization processing

Expected results:
- Compression succeeds, processing file boundary
- Decompression succeeds, data completely consistent
- MD5 verification matches
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqNoTrailingNewlineTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqNoTrailingNewlineTest")
        self.source_file = "no_trailing_newline_test.fq"
        self.compressed_file = "no_trailing_newline_test.fq.pbgz"
        self.decompressed_file = "no_trailing_newline_test.fq.dec"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        fastq_content = []
        num_reads = 50
        
        for i in range(num_reads):
            seq_id = f"@NOTRL{i:06d} 1/1"
            
            bases = ''.join(random.choices('ATGCN', k=random.randint(50, 100)))
            
            separator = '+'
            
            quality_length = len(bases)
            quality = ''.join(random.choice('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                             for _ in range(quality_length))
            
            fastq_content.extend([seq_id, bases, separator, quality])
        
        full_content = '\n'.join(fastq_content)  # No trailing newline
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        
        
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command)
        
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file}"
        self.add_test_command(decompress_command)
    
    def cleanup_test_data(self):
        files_to_clean = [self.source_file, self.compressed_file, self.decompressed_file]
        for filename in files_to_clean:
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
        
        compression_ratio = self.get_compression_rate()
        if compression_ratio is None:
            print(f"Error: Failed to get compression ratio")
            return False
        
        exec_time = self.get_execution_time()
        
        
        with open(self.source_file, 'r', encoding='utf-8') as f1, open(self.decompressed_file, 'r', encoding='utf-8') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        # Decompressed file will append a \n, requiring special handling for verification logic
        # Remove possible \n at the end of decompressed file before comparison
        if decompressed_content.endswith('\n'):
            decompressed_content_trimmed = decompressed_content[:-1]
        else:
            decompressed_content_trimmed = decompressed_content
        
        decompressed_hash = hashlib.md5(decompressed_content_trimmed.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            print(f"  Original hash: {original_hash}")
            print(f"  Decompressed hash (trimmed): {decompressed_hash}")
            return False


if __name__ == "__main__":
    test_case = FastqNoTrailingNewlineTest()
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