"""
testtestusecase：SamWithReferenceTest - SAMfilefilecompresscompressdecompresscompresstesttest（highqualityquantityreferencereferencedatadata）

testtestscenarioscenario：
- testtestpbgztohighqualityquantitySAMfilefileadvanceexecutioncompresscompressanddecompresscompress
- SAMfilefilecontainingcontainpreciseheartsetcountofcomparetodatadata
- useusereferencereferencebasegenegenomerealpresentmostbestcompresscompressvalidresult
- testtestrelatedkeyfieldfield：readnamecall、comparetopositionposition、reflectshootqualityquantity、CIGARfieldcharacterstringequal

useusecommand：
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

referencedatasayclear：
- -rreferencedata：useusereferencereferencebasegenegenomeabqd01.fasta.gz
- testtestdatadatacontainingcontain3itemhighqualityquantityofSAMrecordrecord

Expectedexpectedresultresult：
- compresscompressgeneratesuccess，compresscompressrateveryhigh（geneforreferencereferencematchpair）
- decompresscompressgeneratesuccess，datadatacompletecompleterestorecomplex
- MD5verifyVerifycompletefullonecause
- placehavecanselectfieldfieldcorrectaccurateensurekeep（NM, MD, AS, XS, RGequal）

skilltechniquesayclear：
- testtesthighqualityquantitySAMdatadataofhandlemanageabilityforce
- Verifyverifypbgzhandlemanagecomplexcomplexcanselectfieldfieldofabilityforce
- Expectedexpectedcompresscompressrateveryhigh，genefordatadatacontainingcontainreferencereferencematchpairinformationinformation
"""

"""
testtestusecase：SamWithReferenceTest - SAMhighqualityquantitydatadatacompresscompressdecompresscompresstesttest（referencereferencebasegenegenome）

testtestscenarioscenario：
- testtestpbgztohighqualityquantitySAMfilefileadvanceexecutioncompresscompressanddecompresscompress
- useusereferencereferencebasegenegenomeabqD01.fasta.gzoptimaltransformcompresscompress
- SAMdatadatacontainingcontainpreciseaccurateofreferencereferencematchpair，Expectedexpectedcompresscompressratehigh

useusecommand：
1. ./release-release/bin/pbgz compress {source_file} -o {compressed_file} -r {reference_file}
2. ./release-release/bin/pbgz decompress {compressed_file} -o {decompressed_file} -r {reference_file}

referencedatasayclear：
- -rreferencedata：useusereferencereferencebasegenegenomeoptimaltransformcompresscompress，specialotheristoatcontainingcontainpreciseaccuratematchpairofdatadata
- testtest3itemhighqualityquantityofSAMrecordrecord，positionpositionat1k、7kattachnear

Expectedexpectedresultresult：
- compresscompressgeneratesuccessrateveryhigh（geneforcontainingcontainpreciseaccurateofreferencereferencematchpair）
- decompresscompressgeneratesuccess，datadatacompletebeautyrestorecomplex
- MD5verifyVerifycompletefullonecause
- placehavecanselectfieldfieldcorrectaccurateensurekeep：NM（errorpairdata）、MD（errorpairdetailnode）、AS（comparetodividedata）、XS（timeoptimalcomparetodividedata）、RG（read group）
- Expectedexpectedcompresscompressrateveryhigh，genefordatadatacontainingcontaincompletefullmatchpairofreferencereferenceorderlist

skilltechniquesayclear：
- testtesthighqualityquantitySAMdatadataofhandlemanage
- Verifyverifypbgzhandlemanagecomplexcomplexcanselectfieldfieldofabilityforce
- referencereferencebasegenegenomeabilityextremelargeimproveupgradecompresscompressvalidresult
"""

import os
import sys
import random
import hashlib

# AddparenttargetrecordtoPythonroadpath
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamWithReferenceTest(PBGZTestCase):
    def __init__(self):
        super().__init__("SamWithReferenceTest")
        self.source_file = "with_ref_test.sam"
        self.compressed_file = "with_ref_test.sam.pbgz"
        self.decompressed_file = "with_ref_test.sam.dec"
        self.reference_file = "fa/ABQD01.fasta.gz"
    
    def get_test_files(self) -> tuple:
        """Returnreturntesttestfilefileinformationinformation"""
        return (self.source_file, self.compressed_file, self.decompressed_file)

    def prepare_data(self):
        # GenerategeneratebasebaseofSAMfilefileinsidecapacity
        sam_content = []
        
        # SAMheader
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # Add3itemreadrecordrecord
        num_reads = 3
        for i in range(num_reads):
            # QNAME
            qname = f"SAM_READ_{i:03d}"
            
            # FLAG
            flag = 0
            
            # RNAME
            rname = "ENA|ABQD01000001|ABQD01000001.1"
            
            # POS
            pos = (i + 1) * 1000
            
            # MAPQ
            mapq = 60
            
            # CIGAR
            cigar = "25M"
            
            # RNEXT
            rnext = "*"
            
            # PNEXT
            pnext = 0
            
            # TLEN
            tlen = 0
            
            # SEQ
            seq = ''.join(random.choices('ATGCN', k=25))
            
            # QUAL
            qual = '!' * 25
            
            # Addcanselectfieldfield，withreferencereferencebasegenegenomebasegeneorderlisttesttestonecause
            read_line = f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:25\tAS:i:25\tXS:i:0\tRG:Z:1"
            
            sam_content.append(read_line)
        
        self.original_hash = hashlib.md5('\n'.join(sam_content).encode('utf-8')).hexdigest()
        
        # GenerategenerateSAMfilefile
        with open(self.source_file, 'w') as f:
            for line in sam_content:
                f.write(line + '\n')
        
        
        # Addcompresscompresstesttestcommand
        compress_command = f"./release-release/bin/pbgz compress {self.source_file} -o {self.compressed_file} -r {self.reference_file}"
        self.add_test_command(compress_command, cmd_id=1)
        
        # Adddecompresscompresstesttestcommand
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
        # strictformatCheckcheck：compresscompressanddecompresscompresscommandallmustmustgeneratesuccess
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
            return False
        
        # VerifyverifycompresscompressfilefileiswhetherGenerategenerate
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
        
        # useuseMD5tocompare
        decompressed_hash = hashlib.md5(decompressed_content.encode('utf-8')).hexdigest()
        original_hash = hashlib.md5(original_content.encode('utf-8')).hexdigest()
        
        return original_hash == decompressed_hash


if __name__ == "__main__":
    test_case = SamWithReferenceTest()
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