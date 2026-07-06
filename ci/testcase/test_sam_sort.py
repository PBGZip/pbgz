"""
testtestusecase：SamSortTest - SAMfilefilesortedordersuccessabilitytesttest

testtestscenarioscenario：
- testtestpbgzofsortcommandtoSAMfilefileadvanceexecutionsortedorder
- VerifyverifySAMrecordrecordbyreferencereferenceorderlistpositionpositioncorrectaccuratesortedorder
- testtestsortedorderafterofcompresscompresssuccessability

useusecommand：
1. ./release-release/bin/pbgz sort {input_file} -o {sorted_file}
2. ./release-release/bin/pbgz compress {sorted_file} -o {compressed_file}

referencedatasayclear：
- sortcommand：bySAMheaderheaderofSOfieldfieldadvanceexecutionsortedorder
- -oreferencedata：pointdeterminesortedorderafterofoutputoutputfilefilepositionposition
- SAMfilefilefromunsortedstatestatevariableforcoordinatestatestate

Expectedexpectedresultresult：
- sortcommandgeneratesuccessexecuteexecution，GenerategeneratesortedorderafterofSAMfilefile
- compresscompresscommandgeneratesuccessexecuteexecution
- SAMrecordrecordbypositionpositioncorrectaccuratesortedorder
- sortedorderafterofcompresscompressfilefilequalityquantitygoodgood

skilltechniquesayclear：
- SAMsortedorderiscomparetoafterdivideanalyzeofheavywantstepstep
- needwantcorrectaccuratehandlemanage@HDheaderheaderinofSOfieldfield
- testtestbasebaseofsortsuccessabilityworkworkflowprocess
"""

import os
import sys
import random
import hashlib

# AddparenttargetrecordtoPythonroadpath
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from testcase.pbgz_test_framework import PBGZTestCase


class SamSortTest(PBGZTestCase):
    def __init__(self):

        super().__init__("SamSortTest")
        self.source_file = "sam_unsorted.sam"
        self.sorted_file = "sam_unsorted.sort.sam"
        self.output_dir = "."

    def get_test_files(self) -> tuple:
        """Returnreturntesttestfilefileinformationinformation"""
        return (self.source_file, None, self.sorted_file)

    def prepare_data(self):
        # GenerategenerateonepieceinvalidorderofSAMfilefile
        sam_content = []
        
        # SAMfilefileheaderheader
        sam_content.append("@HD\tVN:1.0\tSO:unsorted")
        sam_content.append("@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000")
        
        # Generategenerateinvalidorderofdatadataexecution（reasonmeaningtypemessypositionposition）
        # firstGenerategeneratehaveorderofdatadata
        ordered_lines = []
        num_reads = 50
        for i in range(num_reads):
            flag = 0
            seq_id = f"SAM_{i:03d}"
            pos = (i + 1) * 100
            mapq = 60
            cigar = "50M"
            rnext = "*"
            pnext = 0
            tlen = 0
            seq = ''.join(random.choices('ATCGN', k=50))
            qual = '!' * 50
            
            line = f"{seq_id}\t{flag}\tENA|ABQD01000001|ABQD01000001.1\t{pos}\t{mapq}\t{cigar}\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\tNM:i:0\tMD:Z:50\tAS:i:50\tXS:i:0\tRG:Z:SORT_TEST"
            ordered_lines.append(line)
        
        # typemessyorderorder
        random.shuffle(ordered_lines)
        sam_content.extend(ordered_lines)
        
        full_content = '\n'.join(sam_content) + '\n'
        self.original_hash = hashlib.md5(full_content.encode('utf-8')).hexdigest()
        
        with open(self.source_file, 'w', encoding='utf-8') as f:
            f.write(full_content)
        
        # Addsortedordertesttestcommand，useuse-oreferencedataclearaccuratepointdetermineoutputoutputfilefilename，-freferencedatastrongcontrolcovercover
        sort_command = f"./release-release/bin/pbgz sort {self.source_file} -o {self.sorted_file} -f"
        self.add_test_command(sort_command, cmd_id=1)
    
    def cleanup_test_data(self):
        files_to_clean = [self.source_file, self.sorted_file]
        for filename in files_to_clean:
            try:
                if os.path.exists(filename):
                    os.remove(filename)
            except Exception as e:
                print(f"Warning: Failed to remove {filename}: {e}")
    
    def verify_expected_results(self) -> bool:
        # Verifyverifysortcommandiswhethergeneratesuccess
        sort_success = False
        for cmd_result in self.test_results.get("commands", []):
            command = cmd_result.get("command", "")
            success = cmd_result.get("success", False)
            if "sort" in command and "pbgz" in command:
                sort_success = success
                
                # likeresultsortlosefail，recordrecorderrorerrorinformationinformationuseatdebugtest
                if not success:
                    error_info = cmd_result.get("error_info", {})
                    print(f"Sort command failed: {error_info}")
        
        if not sort_success:
            return False
            
        if not os.path.exists(self.sorted_file):
            return False
        
        sorted_size = os.path.getsize(self.sorted_file)
        if sorted_size == 0:
            return False
        
        # Verifyverifyfilefileiswhetherrealcorrecthaveorder
        with open(self.sorted_file, 'r') as f:
            lines = f.readlines()
        
        # Checkcheckheaderheaderinformationinformation（headerheaderexecutionuse@openheader）
        header_lines = []
        data_lines = []
        for line in lines:
            if line.startswith('@'):
                header_lines.append(line)
            else:
                data_lines.append(line)
        
        # Verifyverifydatadataexecutioniswhetherbypositionpositionsortedorder
        positions = []
        for line in data_lines:
            if line.strip():
                parts = line.split('\t')
                if len(parts) >= 4:
                    try:
                        pos = int(parts[3])  # number4listisposition
                        positions.append((pos, line.strip()))
                    except (ValueError, IndexError):
                        pass
        
        # Checkcheckpositionpositioniswhetherincreasingincrease
        for i in range(len(positions) - 1):
            if positions[i][0] > positions[i+1][0]:
                return False
        
        # Verifyverifysortedorderafteroffilefileinsidecapacitycompletecompletequality
        with open(self.sorted_file, 'r') as f:
            sorted_content = f.read()
        
        sorted_hash = hashlib.md5(sorted_content.encode('utf-8')).hexdigest()
        
        # byatsortedorderwillchangevariablefilefileinsidecapacityoforderorder，IpluralnotabilitydirectconnectcomparecompareMD5
        # butIpluralcanuseCheckchecksortedorderafteroffilefilecontainingcontainpast tense markercorrectaccurateofheaderheaderanddatadata
        if "@HD\tVN:1.0\tSO:unsorted" not in sorted_content:
            return False
        
        if "@SQ\tSN:ENA|ABQD01000001|ABQD01000001.1\tLN:1000000" not in sorted_content:
            return False
        
        
        return True


if __name__ == "__main__":
    test_case = SamSortTest()
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