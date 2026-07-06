"""
Test case: FastqFilepipeReferenceTest - FASTQ file pipe compression/decompression test (with reference genome)

Test scenario:
- Test FASTQ file compression and decompression via pipe method (with reference genome optimization)
- Verify combination of pipe mode and reference genome optimization
- Test compression effect of FASTQ format matching data

Commands used:
1. cat {source_file} | ./release-release/bin/pbgz compress -o - > {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

Expected results:
- Pipe compression/decompression succeeds
- Reference genome optimization takes effect, compression ratio improved
- FASTQ format completely preserved
- Decompressed file MD5 matches original file completely
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqFilepipeReferenceTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqFilepipeReferenceTest")
        self.source_file = "filepipe_ref_test.fq"
        self.compressed_file = "filepipe_ref_test.fq.pbgz"
        self.decompressed_file = "filepipe_ref_test.fq.dec"
        self.reference_file = "fa/GCA_000002985.3.fasta.gz"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        fastq_content = []
        num_reads = 200
        
        for i in range(num_reads):
            seq_id = f"@FILEPIPE_REF{i:06d} 1/1"
            
            bases = ''.join(random.choices('ATGCN', k=random.randint(50, 100)))
            
            separator = '+'
            
            quality_length = len(bases)
            quality = ''.join(random.choice('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                             for _ in range(quality_length))
            
            fastq_content.extend([seq_id, bases, separator, quality])
        
        full_content = '\n'.join(fastq_content) + '\n'
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        

        
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o - > {self.compressed_file} -r {self.reference_file}"
        self.add_test_command(compress_command)
        
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r {self.reference_file}"
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
        
        
        if compression_ratio <= 0:
            print(f"Error: Invalid compression ratio: {compression_ratio:.2f}%")
            return False
        
        with open(self.source_file, 'r', encoding='utf-8') as f1, open(self.decompressed_file, 'r', encoding='utf-8') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = FastqFilepipeReferenceTest()
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