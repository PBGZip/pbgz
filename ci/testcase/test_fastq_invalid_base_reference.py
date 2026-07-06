"""
Test case: FastqInvalidBaseReferenceTest - FASTQ file invalid base handling test (with reference genome)

Test scenario:
- Test FASTQ files containing invalid bases (non-ATCGN characters) with reference genome optimization
- Verify pbgz's fault tolerance capability and handling strategy for invalid bases
- Test balance between reference genome optimization and data integrity

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

Parameter description:
- FASTQ standard bases are A, T, C, G, N, test illegal characters like Z
- With reference genome, invalid bases may affect compression optimization
- Verify how pbgz handles abnormal data matching with reference genome

Expected results:
- Compression succeeds, handling invalid base data
- Decompression succeeds, data completely preserved (including illegal characters)
- MD5 verification matches completely
"""

import os
import sys
import random
import hashlib

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqInvalidBaseReferenceTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqInvalidBaseReferenceTest")
        self.source_file = "invalid_base_ref_test.fq"
        self.compressed_file = "invalid_base_ref_test.fq.pbgz"
        self.decompressed_file = "invalid_base_ref_test.fq.dec"
        self.reference_file = "fa/GCA_000002985.3.fasta.gz"

    def get_test_files(self) -> tuple:
        """Return test file information"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate FASTQ file containing invalid bases, using reference genome
        fastq_content = []
        num_reads = 80

        for i in range(num_reads):
            # ID line
            seq_id = f"@INV_REF{i:06d} 1/1"

            # Sequence content: mix valid and invalid bases
            valid_bases = ''.join(random.choices('ATGCN', k=random.randint(60, 90)))
            # Insert various invalid characters
            invalid_chars = ['X', 'Z', 'K', 'L', 'M', 'O', 'P', 'Q', 'R', 'S', 'V', 'W', 'Y']
            num_invalid = random.randint(2, 5)
            positions = random.sample(range(1, len(valid_bases)-1), num_invalid)

            seq = list(valid_bases)
            for pos, invalid_char in zip(positions, invalid_chars[:num_invalid]):
                seq[pos] = invalid_char

            seq = ''.join(seq)

            # + line
            separator = '+'

            # Quality score line
            quality_length = len(seq)
            quality = ''.join(random.choice('"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI')
                             for _ in range(quality_length))
            
            fastq_content.extend([seq_id, seq, separator, quality])
        
        full_content = '\n'.join(fastq_content) + '\n'
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        
        
        # Addaddcompresstesttestcommand，Usereferencefile
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file}"
        self.add_test_command(compress_command)
        
        # Addadddecompresstesttestcommand
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r {self.reference_file}"
        self.add_test_command(decompress_command)
    
    def cleanup_test_data(self):
        """clearmanagetesttestusecaseCreateoftemporarytimefile"""
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
        
        
        # tocompareoriginalfileanddecompressfileiswhetheronecause
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
    test_case = FastqInvalidBaseReferenceTest()
    result = test_case.execute()
    test_case.print_results()
    
    if result:
        print("\nTest passed successfully!")
    else:
        print("\nTest failed!")
        # fordebugtestresultresultensuresaveJSONfile
        if not result:
            test_case.save_to_json()
        sys.exit(1)