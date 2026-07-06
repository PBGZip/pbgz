"""
testtestusecase：SamFilepipeTest - SAMfilefilepipelinepipelinecompresscompressdecompresscompresstesttest

testtestscenarioscenario：
- testtestSAMformatformatfilefileThroughthroughUnixpipelinepipelinemethodformatadvanceexecutioncompresscompressdecompresscompresshandlemanage
- VerifyverifySAMformatformatdatadatawithstandardstandardUnixworktoolofconcurrentcapacityquality
- testtestGenerateobjectinformationinformationstudyworkworkflowinSAMdatadataofpipelinepipelinetransmitoutput

useusecommand：
1. cat {source_file} | ./release-release/bin/pbgz compress -o {compressed_file}
2. ./release-release/bin/pbgz decompress {compressed_file} | cat > {decompressed_file}

Expectedexpectedresultresult：
- pipelinepipelinecompresscompressdecompresscompressgeneratesuccess
- SAMformatformatcompletecompleteensurekeep（@header、11listdatadata）
- decompresscompressafterfilefileMD5withoriginalstartfilefilecompletefullonecause
- pipelinepipelinepatternformatsuccessabilitycorrectnormalworkwork
"""

import os
import sys
import random
import hashlib

# AddparenttargetrecordtoPythonroadpath
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamFilepipeTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamFilepipeTest")
        self.source_file = "sam_filepipe_test.sam"
        self.compressed_file = "sam_filepipe_test.sam.pbgz"
        self.decompressed_file = "sam_filepipe_test.sam.dec"

    def get_test_files(self) -> tuple:
        """Returnreturntesttestfilefileinformationinformation"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # GenerategeneratebasebaseofSAMfilefileinsidecapacity
        sam_content = []
        
        # SAMfilefileheaderheader
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # Additemreadrecordrecord
        num_reads = 100
        for i in range(num_reads):
            qname = f"SAM_FILEPIPE_{i:03d}"
            flag = 0
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            pos = (i + 1) * 100
            mapq = 60
            cigar = "30M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=30))
            qual = '!' * 30
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:30\tAS:i:30\tXS:i:0\tRG:Z:FILEPIPE"
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        # Addcompresscompresstesttestcommand：filefileoutputinput，pipelinepipelineoutputoutput
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o - > {self.compressed_file}"
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
        
        
        # tocompareoriginalfilefileanddecompresscompressfilefileiswhetheronecause
        with open(self.source_file, 'r') as f1, open(self.decompressed_file, 'r') as f2:
            original_content = f1.read()
            decompressed_content = f2.read()
        
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        if original_hash == decompressed_hash:
            return True
        else:
            print(f"✗ Filepipe decompression verification failed: Files are different")
            return False


if __name__ == "__main__":
    test_case = SamFilepipeTest()
    result = test_case.execute()
    test_case.print_results()
    
    if result:
        print("\n✓ Filepipe test passed successfully!")
    else:
        print("\n❌ Filepipe test failed!")
        # fordebugtestresultresultensuresaveJSONfilefile
        if not result:
            test_case.save_to_json()
        sys.exit(1)