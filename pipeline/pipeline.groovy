// vim: ts=4:sw=4:expandtab:cindent
/////////////////////////////////////////////////////////////////////////////////
//
// This file is part of Oncopipe.
// 
// Oncopipe is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, under version 3 of the License, subject
// to additional terms compatible with the GNU General Public License version 3,
// specified in the LICENSE file that is part of the Oncopipe distribution.
//
// Oncopipe is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Oncopipe.  If not, see <http://www.gnu.org/licenses/>.
//
/////////////////////////////////////////////////////////////////////////////////
//
// Oncopipe Main Pipeline Script
//
/////////////////////////////////////////////////////////////////////////////////

about title: "Cancer analysis pipeline"

// Load the default configuration
load 'config.groovy'

// Local file can set EXOME_TARGET and ANALYSIS_PROFILES
if(file("../target_regions.txt").exists())  {
    load '../target_regions.txt'
}

requires EXOME_TARGET : """
        The exome target regions. This should be the whole regions targeted
        for capture by the exome kit. Note that the regions for analysis
        may be a subset, but you should always specify the whole exome
        region here.
    """

// All the core pipeline stages in the pipeline
load 'pipeline_stages_config.groovy'

// load wrapper for JAFFA
load 'jaffa_direct_pipeline.groovy'
load 'classifier_stages_config.groovy'

// Legacy from Melbourne Genomics Health Alliance Sequencing Pipeline
sample_metadata_file = correct_sample_metadata_file( args[0] ) // fix syntax issues and update sample_metadata_file

try {
  sample_info = SampleInfo.parse_mg_sample_info(sample_metadata_file)
}
catch (RuntimeException e) {
  sample_info = SampleInfo.parse_sample_info(sample_metadata_file)
}

// We are specifying that each analysis takes place inside a fixed file structure
// where the parent directory is named according to the batch name. Thus we
// can infer the batch name from the name of the parent directory.
// 
// Note: this variable can be overridden by passing a parameter to bpipe in case
// you are running in a different location.
batch = new File("..").canonicalFile.name

// Extract the analysis profiles from the sample information
ANALYSIS_PROFILES = sample_info*.value*.target as Set

samples = sample_info.keySet()

run {
    // Check the basic sample information first
    check_sample_info + // check that fastq files are present
    check_tools + // The JAFFA equivalent is run_check
    call_jaffa_run_check +

    generate_pipeline_id + // make a new pipeline run ID file if required

    // For each analysis profile we run the main pipeline in parallel
    ANALYSIS_PROFILES * [

        set_target_info +  // eg (B_ALL)

        init_analysis_profile +

        // The first phase is to perform transcriptome alignment and fusion
        // detection (modified from JAFFA), for each sample
        // The second parallel phase is to perform genome alignment (using STAR)
        // call the classifier (if RData exists) and variant calling, for each sample
        samples * [
            // input samples 
            set_sample_info +

            // Run the two pipelines in parallel
            [
                // JAFFA pipeline: align to the transcriptome
                // Loaded from jaffa_direct_pipeline.groovy
                // Requires JAFFA 1.09_dev or later
                ~"(.*)_R[1-2][_.].*f(ast)?q.gz" * [ 
                    call_jaffa_direct
                ]
            ,
                // Align to the genome.
                // Call the classifier pipeline
                ~"(.*)_R[1-2][_.].*f(ast)?q.gz" * [
                    fastqc
                    ,
                    classifier_pipeline
                ]
            ]
        ] +
        // TODO: modify to use "AllSorts" (Currently calls
        // first iteration of Classifier R code)
        count_reads_RNA +
        predict_class +
        // TODO: This should really be changed to be per sample
        // Discuss need for all results in one spreadsheet.
        direct_compile_all_results
    ]
}
