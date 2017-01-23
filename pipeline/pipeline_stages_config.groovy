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

set_target_info = {

    doc "Validate and set information about the target region to be processed"

    output.dir="../design"

    branch.batch = batch
    branch.target_name = branch.name
    branch.target_samples = sample_info.grep { it.value.target == target_name }*.value*.sample
    // needed for variant calling and filtering
    branch.target_gene_file = "${output.dir}/${target_name}.genes.txt"

    // needed for classifier
    branch.target_RFmodel_file = "${output.dir}/${target_name}.model.RData"
    branch.target_FPKM_file = "${output.dir}/${target_name}.fpkm.txt"
    branch.target_topgenenames_file = "${output.dir}/${target_name}.topgenenames.RData"

    println "Checking for target random forest file: $target_RFmodel_file"
    produce(target_RFmodel_file) {
        exec """
            if [ -e $BASE/designs/$target_name/${target_name}.model.RData ];
            then
                cp $BASE/designs/$target_name/${target_name}.model.RData $target_RFmodel_file;
            else
                touch $target_RFmodel_file;
            fi;
        """
    }

    println "Checking for target FPKM file: $target_FPKM_file"
    produce (target_FPKM_file) {
        exec """
            if [ -e $BASE/designs/$target_name/${target_name}.fpkm.txt ];
            then
                cp $BASE/designs/$target_name/${target_name}.fpkm.txt $target_FPKM_file;
            else
                touch $target_FPKM_file;
            fi;
        """
    }

    println "Checking for topGeneNames RData file: $target_topgenenames_file"
    produce (target_topgenenames_file) {
        exec """
            if [ -e $BASE/designs/$target_name/${target_name}.topgenenames.RData ];
            then
                cp $BASE/designs/$target_name/${target_name}.topgenenames.RData $target_topgenenames_file;
            else
                touch $target_topgenenames_file;
            fi;
        """
    }

    println "Checking for target gene file: $target_gene_file"
    produce(target_gene_file) {
        exec """
            if [ -e $BASE/designs/$target_name/${target_name}.genes.txt ];
            then
                cp $BASE/designs/$target_name/${target_name}.genes.txt $target_gene_file;
            else
                touch $target_gene_file;
            fi;
        """
    }
}

init_analysis_profile = {
  // This stage is a placeholder to allow individual analysis profiles
  // to perform initialization steps
}

set_sample_info = {

    doc "Validate and set information about the sample to be processed"

    branch.sample = branch.name
    if(sample_info[sample].target != target_name) {
        // This is expected because every file is processed for every target/flagship
        succeed "No files to process for sample $sample, target $target_name"
    }

    // Patient specific variants are not supported yet
    // If they are provided, we should not process the patient at all
//    check {
//        if(sample_info[sample].variantsFile?.trim()) 
//                exec "false" // force failure
//    } otherwise { 
//        succeed """
//             Study $sample is configured with a sample specific variant file. The pipeline currently does 
//             not support sample specific variants. Please remove the variant file from the configuration
//             to allow processing.
//         """.trim().stripIndent() to channel: oncopipe_operator, subject: "Invalid configuration for Study $sample"
//    }
//
    def files = sample_info[sample].files.fastq

    println "Processing input files ${files}"
    forward files
}

check_tools = {
    doc """
        Checks for presence of optional tools and sets appropriate pipeline variables
        to enable or disable corresponding pipeline features
    """

    var UPDATE_FUSION_DB : FUSION_DB,
        ANNOTATION_FUSION_DB : FUSION_DB

    produce("revision.txt") {
        exec """
            git describe --always > $output.txt || true
        """
    }

    if(file(GROOVY_NGS).name in ["1.0.1","1.0"])
        fail "This version of Oncopipe requires GROOVY_NGS >= 1.0.2. Please edit config.groovy to set the latest version of tools/groovy-ngs-utils"

    branch.UPDATE_FUSION_DB = UPDATE_FUSION_DB
    branch.ANNOTATION_FUSION_DB = ANNOTATION_FUSION_DB
}

check_sample_info = {

    doc "Validate basic sample information is correct"

    def missingSummary = []
    for(sample in samples) {

        // Check that FASTQ files begin with the sample name followed by underscore
        def files = sample_info[sample].files.fastq
        if(files.any { !file(it).name.startsWith(sample+"_")}) {
            files.each { println "FASTQ: $it | sample=$sample" }
            fail report('templates/invalid_input.html') to channel: oncopipe_operator, 
                subject: "FASTQ files for sample $sample have invalid file name format", 
                message: "Files $files do not start with the sample name $sample" 
        }

        // Check that all the files specified for the sample exist
        def missingFiles = files.grep { !file(it).exists() }
        if(missingFiles) 
            missingSummary << """
                The following files specified for sample $sample could not be found:\n\n${missingFiles*.center(120).join('\n')}

                Please check that the files in your sample file really exist in the data directory.
            """.stripIndent()

        // Check that file names contain the lane information
        def missingLanes = files.grep { !(it ==~ ".*_L[0-9]*_.*") }
        if(missingLanes) 
            missingSummary << """
                The following files specified for sample $sample do not contain lane information:\n\n${missingLanes*.center(120).join('\n')}

                FASTQ file names are required to contain lane identifiers such as L001, L1 or similar. 
                Please check your input FASTQ and rename it if necessary.
            """

        // Check that file names contain the read number information
        def missingRP = files.grep { !(it ==~ ".*_R[1-2][_.].*f(astq)?.gz\$") }
        if(missingRP) 
            missingSummary << """
                The following files for sample $sample do not contain the read number in the expected format:\n\n${missingRP*.center(120).join('\n')}

                FASTQ file names are required to contain the number of the read from the read pair (1 or 2) 
                in the form '_R1_' or '_R1.'. Please check your input FASTQ and rename it if necessary.
            """
    }

    if(missingSummary) {
        fail missingSummary.join("\n" + ("-" * 120) + "\n")
    }
}


fastqc = {
    doc "Run FASTQC to generate QC metrics for raw reads"
    output.dir = "fastqc"
    transform('.fastq.gz')  to('_fastqc.zip') {
        exec "$FASTQC/fastqc -o ${output.dir} $inputs.gz"
    }
}

check_fastqc = {

    doc "Search for any failures in FastQC output and abort further processing if they are found"

    check {
       // NOTE: we remove per-base-sequence content and
       // per-base-gc-content from examination because Nextera
       // appears to contain natural biases that flag QC failures 
       // here.
       exec """
           cat fastqc/"${sample}"_*fastqc/summary.txt |
               grep -v "Per base sequence content" |
               grep -v "Per base GC content" |
               grep -q 'FAIL' && exit 1

           exit 0
       ""","local"
    } otherwise {
        if(CHECK_FASTQC_FAILURES) {
            succeed report('templates/fastqc_failure.html') to channel: oncopipe_operator, 
                subject: "Sample $sample has failed FastQC Check", 
                file: input.zip
        } else {
            send report('templates/fastqc_failure.html') to channel: oncopipe_operator, 
                subject: "Sample $sample has failed FastQC Check", 
                file: input.zip
        }
    }

    check("FASTQ Format") {
        exec """
            awk -F'\\t' '/Illumina/ { where=match(\$2, /[0-9.]+/); { result=substr(\$2, where, RLENGTH); exit(result<1.7); } }' fastqc/${sample}_*_fastqc/fastqc_data.txt
        ""","local"
    } otherwise {
        println "=" * 100
        println "Sample $sample is encoded using a quality encoding incompatible with this pipeline."
        println "Please convert the data first using maq ill2sanger."
        println "=" * 100

        succeed report('templates/fastqc_failure.html') to channel: oncopipe_operator, 
            subject: "Sample $sample is encoded with incompatible quality scores (Illumina < 1.7)", 
            file: input.zip
    }
}

trim_fastq = {
   output.dir="align"
   if(ADAPTERS_FASTA) {
       filter("trim","trim") {
           exec """
               $TRIMMOMATIC/trimmomatic PE -phred33 
               $input1.gz $input2.gz 
               $output1.gz ${output1.prefix}.unpaired.gz 
               $output2.gz ${output2.prefix}.unpaired.gz 
               ILLUMINACLIP:$ADAPTERS_FASTA:2:40:15 LEADING:3 TRAILING:6 SLIDINGWINDOW:4:15 MINLEN:36
           """
       }
   }
}

cleanup_trim_fastq = {
    if(ADAPTERS_FASTA)
        cleanup "*.fastq.trim.gz"
}

cleanup_initial_bams = {
    cleanup("*.merge.bam", ~".*L00[0-9].bam")
}

cleanup_intermediate_bams = {
    cleanup("*.dedup.bam", "*.realign.bam")
}

// remove spaces from gene lists and point to a new sample metadata file
// note that this isn't run through bpipe
correct_sample_metadata_file = {
    def target = new File( 'results' )
    if( !target.exists() ) {
        target.mkdirs()
    }
    [ "sh", "-c", "python $SCRIPTS/correct_sample_metadata_file.py < $it > results/samples.corrected" ].execute().waitFor()
    return "results/samples.corrected"
}

generate_pipeline_id = {
    doc "Generate a pipeline run ID for this batch"
    output.dir="results"
    produce("run_id") {
      exec """
        python $SCRIPTS/update_pipeline_run_id.py --id $ID_FILE --increment True > $output
      """
    }
   // This line is necessary on some distributed file systems (e.g. MCRI) to ensure that
   // files get synced between nodes
   file("results").listFiles()
   run_id = new File('results/run_id').text.trim()
}

create_sample_metadata = {
    doc "Create a new samples.txt file that includes the pipeline ID"
    requires sample_metadata_file : "File describing meta data for pipeline run (usually, samples.txt)"

    output.dir="results"
    produce("results/samples.meta") {
      exec """
          python $SCRIPTS/update_pipeline_run_id.py --id results/run_id --parse True < $sample_metadata_file > results/samples.meta
      """
    }
}
