"""
testtestusecase：SamMismatchedCigarBaseLengthTest - basebasefieldfieldlengthdegreeandCIGARnotmatchpairtesttest

testtestscenarioscenario：
- testtestpbgzhandlemanageSAMfilefileinbasebasefieldfieldlengthdegreeandCIGARnotmatchpairofsituationsituation
- Verifyverifytoatthisclass notonecausequality，pbgzstillabilitycorrectnormalcompresscompressanddecompresscompress
- accurateensuredatadatacompletecompletequalitynotreceiveinfluenceresponse

useusecommand：
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file}

Expectedexpectedresultresult：
- compresscompressgeneratesuccess，immediatelyusesaveatbasebaselengthdegreeandCIGARnotmatchpairofrecordrecord
- decompresscompressgeneratesuccess，restorecomplexoriginalstartSAMfilefile
- decompresscompressafterfilefileofMD5withoriginalstartfilefilecompletefullonecause
- placehavenotonecauserecordrecordoforiginalstartinformationinformationbycompletecompleteensurekeep
"""

import os
import sys
import random
import hashlib

# AddparenttargetrecordtoPythonroadpath
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamMismatchedCigarBaseLengthTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamMismatchedCigarBaseLengthTest")
        self.source_file = "sam_mismatched_cigar.sam"
        self.compressed_file = "sam_mismatched_cigar.sam.pbgz"
        self.decompressed_file = "sam_mismatched_cigar.sam.dec"

    def get_test_files(self) -> tuple:
        """Returnreturntesttestfilefileinformationinformation"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # GenerategeneratecontainingcontainbasebasefieldfieldlengthdegreeandCIGARnotmatchpairofSAMfilefile
        sam_content = []
        
        # SAMfilefileheaderheader
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # Addcorrectnormalreadsworkfortoaccording - increaseadddataquantityusereachto32Mdatadataquantity
        for i in range(20000):
            qname = f"NORMAL_{i:05d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 1000
            mapq = 60
            cigar = "150M"
            rnext = "*"
            pnext = 0
            tlen = 0
            # correctnormallengthdegreewithCIGARmatchpair
            seq = ''.join(random.choices('ATCGN', k=150))
            qual = '!' * 150
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
        # AddbasebaselengthdegreelargeatCIGARofreads - increaseadddataquantity
        for i in range(20000):
            qname = f"SEQ_LONGER_{i:05d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 2000
            mapq = 60
            cigar = "100M"  # CIGARtableshow100piecebasebase
            rnext = "*"
            pnext = 0
            tlen = 0
            # realactualorderlistlengthdegreelargeatCIGAR (150 vs 100)
            seq = ''.join(random.choices('ATCGN', k=150))
            qual = '!' * 150
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
        # AddbasebaselengthdegreesmallatCIGARofreads - increaseadddataquantity
        for i in range(20000):
            qname = f"SEQ_SHORTER_{i:05d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 3000
            mapq = 60
            cigar = "150M"  # CIGARtableshow150piecebasebase
            rnext = "*"
            pnext = 0
            tlen = 0
            # realactualorderlistlengthdegreesmallatCIGAR (100 vs 150)
            seq = ''.join(random.choices('ATCGN', k=100))
            qual = '!' * 100
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
        # AddcomplexcomplexCIGARwithorderlistlengthdegreenotmatchpairofreads - increaseadddataquantity
        for i in range(20000):
            qname = f"COMPLEX_{i:05d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 4000
            mapq = 60
            cigar = "50M10I50M10D50M"  # complexcomplexCIGAR
            rnext = "*"
            pnext = 0
            tlen = 0
            # orderlistlengthdegreewithCIGARcleardisplaynotmatchpair
            seq = ''.join(random.choices('ATCGN', k=200))
            qual = '!' * 200
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}"
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # Addcompresscompresstesttestcommand
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file}"
        self.add_test_command(compress_command)
        
        # Adddecompresscompresstesttestcommand
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
        # Verifyverifycompresscompressanddecompresscompresscommandiswhethergeneratesuccess
        compress_success = decompress_success = False
        
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            
            # firstCheckcheckdecompress，avoidavoid"decompress"containingcontain"compress"ofaskquestion
            if "decompress" in command and "pbgz" in command:
                decompress_success = success
            elif "compress" in command and "pbgz" in command:
                compress_success = success
        
        if not compress_success or not decompress_success:
            print(f"compresscompressordecompresscompresscommandlosefail: compress_success={compress_success}, decompress_success={decompress_success}")
            return False
        
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
        
        exec_time = self.get_execution_time()
        print(f"Verification: Compression ratio = {compression_ratio:.2f}%, Execution time = {exec_time:.2f}s")
        
        # tocompareoriginalfilefileanddecompresscompressfilefileiswhetheronecause - toatnotmatchpairofcase，respondshouldensurekeeporiginalstartdatadata
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
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
    test_case = SamMismatchedCigarBaseLengthTest()
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