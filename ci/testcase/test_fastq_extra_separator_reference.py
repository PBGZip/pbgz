"""
Test case: FastqExtraSeparatorReferenceTest - FASTQ file extra separator handling test (with reference genome)

Test scenario:
- Test pbgz handling FASTQ files containing extra separators (with reference genome optimization)
- Verify fault tolerance of non-standard formats combined with reference genome usage
- Test data integrity in extreme cases

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

Parameter description:
- FASTQ format requires strict separators, test fault tolerance for extra separators
- -r parameter: Use reference genome for optimized compression
- Verify combination of file format fault tolerance and reference genome optimization

Expected results:
- Compression succeeds, optimized via reference genome
- Decompression succeeds, data completely restored
- MD5 verification matches completely
"""

import os
import sys
import random
import hashlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqExtraSeparatorReferenceTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqExtraSeparatorReferenceTest")
        self.source_file = "extra_separator_ref.fq"
        self.compressed_file = "extra_separator_ref.fq.pbgz"
        self.decompressed_file = "extra_separator_ref.fq.dec"
        self.reference_file = "fa/GCA_000002985.3.fasta.gz"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        fastq_content = []
        num_reads = 200
        
        for i in range(num_reads):
            seq_id = f"@TEST{i:06d} 1/1"
            
            bases = ''.join(random.choices('ATGCN', k=random.randint(50, 100)))
            
            if i == 50:
                separator = '+'
                fastq_content.extend([seq_id, bases, separator])
                fastq_content.append('+extra')
                quality_length = len(bases)
                quality = ''.join(random.choice('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                                 for _ in range(quality_length))
                fastq_content.append(quality)
            else:
                separator = '+'
                fastq_content.extend([seq_id, bases, separator])
                quality_length = len(bases)
                quality = ''.join(random.choice('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                                 for _ in range(quality_length))
                fastq_content.append(quality)
        
        full_content = '\n'.join(fastq_content) + '\n'
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        

        
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file}"
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
    test_case = FastqExtraSeparatorReferenceTest()
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