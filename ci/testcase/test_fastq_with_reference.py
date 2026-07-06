"""
testtestusecase：FastqWithReferenceTest - FASTQfilefilecompresscompressdecompresscompresstesttest（withreferencereferencebasegenegenome）

testtestscenarioscenario：
- testtestuseusereferencereferencebasegenegenomeofFASTQfilefilecompresscompress
- VerifyverifyreferencereferencebasegenegenometoFASTQcompresscompressofoptimaltransformvalidresult
- testtestFASTQformatformatwithreferencereferencebasegenegenomeofconcurrentcapacityquality

useusecommand：
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

Expectedexpectedresultresult：
- compresscompressgeneratesuccess，benefitusereferencereferencebasegenegenomeoptimaltransformcompresscompresscompare
- decompresscompressgeneratesuccess，datadatacompletebeautyrestorecomplex
- MD5verifyVerifycompletefullonecause
"""

import os
import sys
import random
import hashlib

# AddparenttargetrecordtoPythonroadpath
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class FastqWithReferenceTest(PBGZTestCase):
    def __init__(self):

        super().__init__("FastqWithReferenceTest")
        self.source_file = "with_reference_test.fq"
        self.compressed_file = "with_reference_test.fq.pbgz"
        self.decompressed_file = "with_reference_test.fq.dec"
        self.reference_file = "fa/GCA_000002985.3.fasta.gz"  # useusereferencereferencebasegenefilefile

    def get_test_files(self) -> tuple:
        """Returnreturntesttestfilefileinformationinformation"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # GenerategenerateuseusereferencereferencebasegeneofFASTQfilefiledatadata
        fastq_content = []
        num_reads = 200  # Generategenerate200itemorderlist
        
        for i in range(num_reads):
            # IDexecution：@openheaderoforderliststandardrecognizecharacter
            seq_id = f"@TEST{i:06d} 1/1"
            
            # orderlistinsidecapacity：onlycontainingcontainATGCNoffollowmachineorderlist
            bases = ''.join(random.choices('ATGCN', k=random.randint(50, 100)))
            
            # + execution：divideseparatecharacter
            separator = '+'
            
            # qualityquantitydividedataexecution：Generategeneratewithorderlistlengthdegreematchpairofqualityquantitydividedata（useusenormalseeofPhredqualityquantitydividedatamodelsurround）
            quality_length = len(bases)
            # useuseASCII 33-73modelsurroundinsideoffieldcharacter（torespondPhredqualityquantitydividedata0-40）
            quality = ''.join(random.choice('\"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHI') 
                             for _ in range(quality_length))
            
            fastq_content.extend([seq_id, bases, separator, quality])
        
        full_content = 'n'.join(fastq_content) + 'n'  # correctnormalFASTQfilefileuseexchangeexecutioncharacterresulttail
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        
        
        # Addcompresscompresstesttestcommand，useusereferencereferencebasegenefilefile
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file}"
        self.add_test_command(compress_command, cmd_id=1)
        
        # Adddecompresscompresstesttestcommand，samelikeuseusereferencereferencebasegenefilefile
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r {self.reference_file}"
        self.add_test_command(decompress_command, cmd_id=2)
    
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
        
        
        # Verifyverifycompresscompresscompareiswhetherhavevalid
        if compression_ratio <= 0:
            print(f"Error: Invalid compression ratio: {compression_ratio:.2f}%")
            return False
        
        # tocompareoriginalfilefileanddecompresscompressfilefileiswhetheronecause
        with open(self.source_file, 'r', encoding='utf-8') as f1, open(self.decompressed_file, 'r', encoding='utf-8') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        # useusehashhopevaluetocompare
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = FastqWithReferenceTest()
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