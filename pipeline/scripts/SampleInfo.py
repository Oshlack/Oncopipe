import csv
import os
import re
import sys
import argparse

"""
# BEGIN SAMPLE_INFO TESTS
    >>> si = SampleInfo()
    >>> si.parse_sample_info('samples.txt')
    >>> si = extract_sample_info_from_filename('L10248_S12_L001_R2_001.fastq.gz')
    >>> si = extract_sample_info_from_filename('GHS025_D1EH2ACXX_ATCACG_L001_R1.fastq.bam')
    >>> si = extract_sample_info_from_filename('008_11_H14NGADXX_L1.1.fastq.trim.atrim.reorder.bam')
    >>> si = extract_sample_info_from_filename('MLM-11_ML150083_ML150089_150626_H05ELAFXX_NS500543_L004_R2.fastq.gz')

# END SAMPLE_INFO TESTS
"""

FIELD_TYPES = [ ('Batch', str),
                ('Sample_ID', str),
                ('DNA_Tube_ID', str),
                ('Sex', str),
                ('DNA_Concentration', str),
                ('DNA_Volume', str),
                ('DNA_Quantity', str),
                ('DNA_Quality', str),
                ('DNA_Date', str),
                ('Cohort', str),
                ('Sample_Type', str),
                ('Fastq_Files', str),
                ('Prioritised_Genes', str),
                ('Consanguinity', str),
                ('Variants_File', str),
                ('Pedigree_File', str),
                ('Ethnicity', str),
                ('VariantCall_Group', str),
                ('Capture_Date', str),
                ('Sequencing_Date', str),
                ('Mean_Coverage', str),
                ('Duplicate_Percentage', str),
                ('Machine_ID', str),
                ('DNA_Extraction_Lab', str),
                ('Sequencing_Lab', str),
                ('Exome_Capture', str),
                ('Library_Preparation', str),
                ('Barcode_Pool_Size', str),
                ('Read_Type', str),
                ('Machine_Type', str),
                ('Sequencing_Chemistry', str),
                ('Sequencing_Software', str),
                ('Demultiplex_Software', str),
                ('Hospital_Centre', str),
                ('Sequencing_Contact', str),
                ('Pipeline_Contact', str),
                ('Notes', str) ]

# BEGIN  SampleInfo
class SampleInfo:
    """Sample info class """
    sample_id = ""
    batch = ""
    cohort = ""
    unit = ""
    lane = ""
    machine_id = ""
    sequencing_date = ""
    fastq_files = []
    filename = ""

    '''constructor with dictionary of fields'''
    def __init__(self, field_data={}):
        self.fastq_files = field_data['Fastq_Files'].split(',')
        self.sample_id = field_data['Sample_ID']
        self.batch = field_data['Batch']
        self.cohort = field_data['Cohort']

    def __repr__(self):
        #return "[ sample:'%s', files:%s, target:'%s' ]" % (self.sample_id, self.fastq_files, self.cohort)
        return "'%s': %s" % (self.sample_id, self)
        #return "%s:[sample:'%s', files:[fastq:%s], target:'%s' ]" % (self.sample_id, self.sample_id, self.fastq_files, self.cohort)

    def __str__(self):
        return "[sample:'%s', files:[fastq:%s], target:'%s']" % (self.sample_id, self.fastq_files, self.cohort)

# END SampleInfo

''' parses the sample info txt file '''
def parse_sample_info(filename):
    si_list = []
    with open(filename) as f:
        tdf = csv.DictReader(f, delimiter='\t')
        row = {}
        for row in tdf:
            row.update((key, conversion(row[key]))
                           for key, conversion in FIELD_TYPES)
            si_list.append(SampleInfo(row))
    return si_list

''' extracts sample information from the filename '''
def extract_sample_info_from_filename(filename):
    basename = os.path.basename(filename)
    si = SampleInfo()
    unit, sample, lane = "","",""

    '''remove the index bar code if it is in the filename'''
    result = re.sub(r"_([ATCG]{6,8})_","",basename)

    '''if there is a run identifier, remove it'''
    result = re.sub(r"_RUN[0-9]*_","_",basename)

    laneMatch = re.search(r"_(L[0-9]{1,3})[._]",result)
    if (laneMatch):
        si.lane = laneMatch.group(1)

    '''Illumina machine names end with XX'''
    machineNameMatch = re.search(r"_([A-Z0-9]*XX)_",result)

    if (machineNameMatch):
        si.machine_id =  machineNameMatch.group(1)
        si.unit = si.lane + "." + machineNameMatch.group(1)
        si.sample_id = result[0:machineNameMatch.start(0)]
    else:
        if(laneMatch):
            si.sample_id = result[0:laneMatch.start(0)]
        else:
            si.sample_id = result[0:result.find("_")]

    return si

#for si in parse_sample_info('samples.txt'):
#    print(si)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse sample info')
    parser.add_argument('input_file', help='sample input file (tab-delimited)')
    args = parser.parse_args()
    si_list = parse_sample_info(args.input_file)
    sys.stdout.write("SAM_INFO = %s" % si_list)
