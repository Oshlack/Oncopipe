#!/usr/bin/env python

'''
/////////////////////////////////////////////////////////////////////////////////
//
// This file is part of Cpipe.
// 
// Cpipe is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, under version 3 of the License, subject
// to additional terms compatible with the GNU General Public License version 3,
// specified in the LICENSE file that is part of the Cpipe distribution.
//
// Cpipe is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Cpipe.  If not, see <http://www.gnu.org/licenses/>.
//
/////////////////////////////////////////////////////////////////////////////////
'''
import argparse
import datetime
import glob
import os
import os.path
import sys
__version_info__ = (0,0,1)
__version__ = '.'join(map(str, __version_info__))
__author__ = "Rebecca Louise Evans"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Manage sample info')
    parser.add_argument('data_files', help='data files in fastq format')
    parser.add_argument('--batch_id', required=True, help='batch to which samples belong')
    parser.add_argument('--disease', required=True, help='disease cohort to which samples belong')
    parser.add_argument('--simple', required=False, help='use simple file format')
    args = parser.parse_args()
    sample_info_re = re.compile(r"^[ \t]+$", re.M)

#TODO: what are the requirements?
class SampleInfo(object):
    batch_id = None
    sample_id = None
    fastq_files = None
    metadata = None
    library = None
    unit = None
    lane = None
    machineName = None

    def __init__(self, batch_id, metadata):
        self.batch_id = {}
        self.metadata = {}

def parse_sample_name

CliBuilder cli = new CliBuilder()
cli.with {
    batch "batch to which samples belong", args:1, required: true
    disease "disease cohort to which samples belong", args:1, required: true
    simple "use simple file format", args:0, required: false
}

opts = cli.parse(args)
if(!opts) {
  System.exit(0)
}

samples = SampleInfo.fromFiles(opts.arguments() as List)
samples.each { it.value.batch = opts.batch; it.value.target = opts.disease }

if (opts.simple) {
  println(
    samples*.value*.toTsv().join("\n")
  )
}
else {
  // print header
  println( [ "Batch", "Sample_ID", "DNA_Tube_ID", "Sex", "DNA_Concentration", "DNA_Volume", "DNA_Quantity", "DNA_Quality", "DNA_Date", "Cohort", "Sample_Type", "Fastq_Files", "Prioritised_Genes", "Consanguinity", "Variants_File", "Pedigree_File", "Ethnicity", "VariantCall_Group", "Capture_Date", "Sequencing_Date", "Mean_Coverage", "Duplicate_Percentage", "Machine_ID", "DNA_Extraction_Lab", "Sequencing_Lab", "Exome_Capture", "Library_Preparation", "Barcode_Pool_Size", "Read_Type", "Machine_Type", "Sequencing_Chemistry", "Sequencing_Software", "Demultiplex_Software", "Hospital_Centre", "Sequencing_Contact", "Pipeline_Contact", "Notes" ].join( '\t' ) )
  for (sample in samples) {
    fastq = sample.value.files.collect { it.key == "all" ? [] : it.value }.flatten().join(",")
    geneCategories = sample.value.geneCategories.collect { it.key + ":" + it.value.join(",") }.join(" ")
    println( [ sample.value.batch, sample.value.sample, "", sample.value.sex.encode(), "", "", "", "", "", sample.value.target, "Normal", fastq, geneCategories, "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "", "" ].join( '\t' ) )
  }
}
