"""
Test case: SamFileTest - Basic SAM file compression/decompression test (with reference genome)

Test scenario:
- Test pbgz's compression and decompression for SAM (Sequence Alignment Map) files
- SAM file contains high-quality genomic alignment data
- Use reference genome for compression to improve compression efficiency

Commands used:
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

Parameter description:
- -r parameter: Specify reference genome file for optimizing compression ratio
- Test data contains 1000 SAM records, positions ranging from 97-900

Expected results:
- Compression command succeeds, generating .pbgz file
- Decompression command succeeds, restoring original SAM file
- Decompressed file MD5 matches original file completely
- compresscompressrategoodgood（geneforuseusereferencebasegenegenome）
- SAMheaderheaderanddatadataexecutioncompletecompleteensurekeep

skilltechniquesayclear：
- SAMformatformatisbasegenegenomecomparetoofstandardstandardformatformat
- indexusereferencebasegenegenomeabilitylargelevelimproveupgradecompresscompressvalidresult
- Verifyverifycomparetoinformationinformation、qualityquantityvalue、canselectfieldfieldequalofcompletecompletequality
"""

import os
import sys
import random
import hashlib

# Add parent directory to Python path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamFileTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamFileTest")
        self.source_file = "sam_only_test.sam"
        self.compressed_file = "sam_only_test.sam.pbgz"
        self.decompressed_file = "sam_only_test.sam.dec"

    def get_test_files(self) -> tuple:
        """Returnreturntesttestfilefileinformationinformation"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # Generate SAM file data
        sam_content = []

        # SAM file header
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000")
        
        # Addlargequantityreadgetrecordrecordadvanceexecutioncompresscompresstesttest
        num_reads = 1000
        for i in range(num_reads):
            # QNAME
            qname = f"SAM_READ_{i:04d}"
            
            # FLAG
            flag = 0
            
            # RNAME
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            
            # POS
            pos = (i % 900) + 1
            
            # MAPQ
            mapq = 60
            
            # CIGAR
            cigar = "50M"
            
            # RNEXT
            rnext = "*"
            
            # PNEXT
            pnext = 0
            
            # TLEN
            tlen = 0
            
            # SEQ
            seq = ''.join(random.choices('ATCGN', k=50))
            
            # QUAL
            qual = '!' * 50
            
            # canselectfieldfield
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:50\tAS:i:50\tXS:i:0\tRG:Z:1"
            
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        # Generate SAM file
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # Addcompresscompresstesttestcommand
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -r fa/ABQD01.fasta.gz"
        self.add_test_command(compress_command, cmd_id=1)
        
        # Adddecompresscompresstesttestcommand
        decompress_command = f"./release-release/bin/pbgz decompress {self.compressed_file} -o {self.decompressed_file} -r fa/ABQD01.fasta.gz"
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
        # strictformatCheckcheck：compresscompressanddecompresscompresscommandallmustmustgeneratesuccess
        compress_success = decompress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
                    # Check decompress first to avoid "decompress" containing "compress" issue
            if "decompress" in command and "pbgz" in command:
                decompress_success = success
            elif "compress" in command and "pbgz" in command:
                compress_success = success
        
        if not compress_success or not decompress_success:
            print(f"Compression or decompression command failed: compress_success={compress_success}, decompress_success={decompress_success}")
            return False

        # Verify compressed file is generated
        if not os.path.exists(self.compressed_file):
            return False
        
        # VerifyverifydecompresscompressfilefileiswhetherGenerategenerate  
        if not os.path.exists(self.decompressed_file):
            return False
        
        # Verifyverifycompresscompresscompareiswhethermatchmanage
        compression_ratio = self.get_compression_rate()
        if compression_ratio is None:
            return False
        
        # tocompareoriginalfilefileanddecompresscompressfilefileiswhetheronecause
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        # useusehashhopevaluetocompare
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        return original_hash == decompressed_hash


if __name__ == "__main__":
    test_case = SamFileTest()
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