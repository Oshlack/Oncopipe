//  This file is a wrapper for JAFFA
// Requires JAFFA 1.09_dev or later
// Ensure codeBase is modified in JAFFA
load "$BASE/tools/JAFFA/JAFFA_stages.groovy"

//fastqInputFormat=~"(.*)_R[0-9][_.].*fastq.gz"
jaffa_output="jaffa/"

// override bowtie2 mapParams to include a metrics output file
//mapParams=$mapParams+" --met-file $input.metrics"

// override prepare_reads to include metrics
direct_prepare_reads = {
    doc "Prepare reads"
    output.dir=jaffa_output+branch
    if (input.size() == 1) {  // single reads
        produce(branch+"_filtered_reads.fastq.gz",
                branch+"_leftover_reads.fastq.gz",
                branch+"_trans.metrics",
                branch+"_masked.metrics") {
            exec """
                $trimmomatic SE -threads $threads -phred$scores $input.gz
                    ${output.dir}/${branch}_trim.fastq
                    LEADING:$minQScore TRAILING:$minQScore MINLEN:$minlen ;
                $bowtie2 $mapParams --very-fast --met-file ${output.dir}/${branch}_trans.metrics
                    --al-gz $output1
                    --un ${output.dir}/temp_trans_unmap_reads.fastq
                    -p $threads -x $transFasta.prefix
                    -U ${output.dir}/${branch}_trim.fastq
                    -S ${output.dir}/${branch}_trans.sam;
                $bowtie2 $mapParams --very-fast --met-file ${output.dir}/${branch}_masked.metrics
                    --un-gz $output2 -p $threads -x $maskedGenome
                    -U ${output.dir}/temp_trans_unmap_reads.fastq 
                    -S ${output.dir}/${branch}_masked.sam;
                cat $output2 >> $output1 ;
                rm ${output.dir}/temp_trans_unmap_reads.fastq ${output.dir}/${branch}_trim.fastq
            ""","direct_prepare_reads"
        }
    } else {  // paired reads
        produce(branch+"_filtered_reads.fastq.1.gz",
                branch+"_filtered_reads.fastq.2.gz",
                branch+"_leftover_reads.fastq.1.gz",
                branch+"_leftover_reads.fastq.2.gz",
                branch+"_trans.metrics",
                branch+"_masked.metrics") {
                // need to check here for whether the files are zipped - FIX
                //trim & fix the file names so Trinity handles the paired-ends reads correctly
            exec """
                $trimmomatic PE -threads $threads -phred$scores $input1 $input2
                    ${output.dir}/tempp1.fq /dev/null
                    ${output.dir}/tempp2.fq /dev/null
                    LEADING:$minQScore TRAILING:$minQScore MINLEN:$minlen;
                function fix_ids {
                    cat \$1 |
                    awk -v app=\$2
                        'BEGIN{ i=0 }{
                        if(i==0) print \$1 \"/\" app ;
                        else print \$1 ;
                        i++ ;
                        if(i==4) i=0 }'
                    2>/dev/null
                ; } ;
                fix_ids ${output.dir}/tempp1.fq 1 > ${output.dir}/${branch}_trim1.fastq ;
                fix_ids ${output.dir}/tempp2.fq 2 > ${output.dir}/${branch}_trim2.fastq ;
                rm ${output.dir}/tempp1.fq ${output.dir}/tempp2.fq ;

                $bowtie2 $mapParams --very-fast --met-file ${output.dir}/${branch}_trans.metrics
                    --al-conc-gz ${output1.prefix.prefix}.gz
                    --un-conc ${output.dir}/temp_trans_unmap_reads.fastq
                    -p $threads -x $transFasta.prefix
                    -1 ${output.dir}/${branch}_trim1.fastq
                    -2 ${output.dir}/${branch}_trim2.fastq
                    -S ${output.dir}/${branch}_trans.sam;
                $bowtie2 $mapParams --very-fast --met-file ${output.dir}/${branch}_masked.metrics
                    --un-conc-gz ${output3.prefix.prefix}.gz
                    -p $threads -x $maskedGenome
                    -1 ${output.dir}/temp_trans_unmap_reads.1.fastq
                    -2 ${output.dir}/temp_trans_unmap_reads.2.fastq
                    -S ${output.dir}/${branch}_masked.sam;
                cat $output3 >> $output1 ;
                cat $output4 >> $output2 ;
                rm  ${output.dir}/temp_trans_unmap_reads.1.fastq
                    ${output.dir}/temp_trans_unmap_reads.2.fastq
                    ${output.dir}/${branch}_trim1.fastq
                    ${output.dir}/${branch}_trim2.fastq;
            ""","direct_prepare_reads"
        }
    }
}

direct_get_unmapped = {
    doc "Get Unmapped"
    output.dir=jaffa_output+branch
    produce(branch+".fasta", branch+"_discordant_pairs.bam") {
        from("*_leftover_reads*.gz") {
            exec """
               $bowtie2 -k1 -p $threads --un ${output.dir}/unmapped.fastq -x $transFasta.prefix
                   -U $input1,$input2 |
               $samtools view -F 4 -S -b - |
               $samtools sort - $output2.prefix ;
               $samtools index $output2 ;
            ""","direct_get_unmapped"
        }
        exec """
            $reformat in=${output.dir}/unmapped.fastq out=${output.dir}/temp.fasta threads=$threads ;
            $dedupe in=${output.dir}/temp.fasta out=$output1 threads=$threads absorbcontainment=f ;
            rm ${output.dir}/temp.fasta ${output.dir}/unmapped.fastq 2> /dev/null
        ""","direct_get_unmapped"
    }
}

direct_align_reads_to_annotation = {
    doc "Align reads to annotation"
    output.dir=jaffa_output+branch
    var CUTOFF_READ_LENGTH : 100
    var DEFAULT_TILE_SIZE : 15
    var DEFAULT_LARGE_TILE_SIZE : 18
    var SAMPLE_SIZE : 1000
    produce(branch+".psl") {
        // find out the average read length so we can set blat's tileSize accordingly
        // we just use the first few thousand from the fasta file as a sample
        from(".fasta") {
            exec """
                SEQ_COUNT=`head -n $SAMPLE_SIZE $input | grep "^>" | wc -l` ;
                SUM_READ_LENGTHS=`head -n $SAMPLE_SIZE $input | grep -v ">" | tr -d "\\n" | wc --chars` ;
                AVERAGE_READ_LENGTH=`expr $SUM_READ_LENGTHS / $SEQ_COUNT` ;
                if [ $readTile -eq "0" ] ; then
                    if [ $AVERAGE_READ_LENGTH -le $CUTOFF_READ_LENGTH ] ; then
                        readTile=$DEFAULT_TILE_SIZE ;
                    else readTile=$DEFAULT_LARGE_TILE_SIZE ;
                    fi ;
                else readTile=$readTile ;
                fi ;
                echo "Using tileSize of \$readTile" ;
                function run_blat {
                    time $blat $transFasta $input -minIdentity=$minIdTrans -minScore=$minScore -tileSize=\$1
                        -maxIntron=$maxIntron $output 2>&1 | tee ${output.dir}/log_blat ;
                } ;
                run_blat \$readTile;
                `### test for the Blat tileSize bug (version 35) ###` ;
                if [[ `cat ${output.dir}/log_blat` == *"Internal error genoFind.c"* ]] ; then
                   echo "Blat error with tileSize=\$readTile" ;
                   echo "Let's try again with tileSize=15" ;
                   run_blat $DEFAULT_TILE_SIZE;
                fi ;
            ""","direct_align_reads_to_annotation"
        }
    }
}

direct_filter_transcripts = {
    doc "Filter transcripts"
    output.dir=jaffa_output+branch
    produce(branch+".txt") {
        from(".psl") {
            exec """
                time $R --vanilla --args $input $output $gapSize $transTable <
                    $R_filter_transcripts_script &> ${output.dir}/log_filter
            ""","direct_filter_transcripts"
        }
    }
}

// verify the grep command here .. -Fx ?? should just be -F
// grep -v "^--" should be "^--$" precisely
// what does reformat do 
// rework the from clause 
// input2 input3 => input1 input2
// from (.psl, .txt, .fasta)
direct_extract_fusion_sequences = {
    doc "Extract fusion sequences"
    output.dir=jaffa_output+branch
    produce(branch+".fusions.fa") {
        from(".txt", ".fasta") {
            exec """
                cat $input1 | awk '{print \$1}' | sed 's/^/>/g' > ${output}.temp ;
                $reformat in=$input2 out=stdout.fasta fastawrap=0 | awk '{print \$1}' |
                grep -F -A1 -f ${output}.temp | grep -v "^\\-\\-" > $output ;
                rm ${output}.temp ;
            ""","direct_extract_fusion_sequences"
        }
    }
}

//input1
direct_align_transcripts_to_genome = {
    doc "Align transcripts to genome"
    output.dir=jaffa_output+branch
    produce(branch+"_genome.psl") {
        from (".fusions.fa") {
            exec """
                $blat $genomeFasta $input -minScore=$minScore $output 2>&1 |
                tee ${output.dir}/log_genome_blat
            ""","direct_align_transcripts_to_genome"
        }
    }
}

direct_make_simple_reads_table = {
    doc "Make simple reads table"
    output.dir=jaffa_output+branch
    produce(branch+".reads") {
        from(".txt", "*_discordant_pairs.bam") {
            exec """
                $samtools view -H $input2 | grep "@SQ" | cut -f2 | sed 's/SN://g' > ${output.dir}/temp_gene_ids ;
                $R --no-save --args $input1 $transTable ${output.dir}/temp_gene_ids ${output.dir}/paired_contigs.temp
                    < $R_get_spanning_reads_direct_script1 ;
                function get_spanning_pairs {
                    gene=`echo \$1 | cut -f1 -d"?"` ;
                    g1=`echo \$1 | cut -f2 -d"?"` ;
                    g2=`echo \$1 | cut -f3 -d "?"` ;
                    $samtools view $input2 \$g1 | cut -d "/" -f1 | sort -u > ${output.dir}/g1 ;
                    $samtools view $input2 \$g2 | cut -d "/" -f1 | sort -u > ${output.dir}/g2 ;
                    left=`cat ${output.dir}/g1 | wc -l` ;
                    right=`cat ${output.dir}/g2 | wc -l` ;
                    both=`cat ${output.dir}/g1 ${output.dir}/g2 | sort -u | wc -l` ;
                    echo -e "\$gene\t\$(( \$left + \$right - \$both ))" ;
                } ;
                cat ${output.dir}/paired_contigs.temp | while read line ; do
                    get_spanning_pairs "\$line" >> ${output.dir}/spanning_pair_counts.temp ;
                done ;
                $R --no-save --args $input1 ${output.dir}/spanning_pair_counts.temp $output < $R_get_spanning_reads_direct_script2 ;
                rm ${output.dir}/temp_gene_ids ${output.dir}/spanning_pair_counts.temp ${output.dir}/paired_contigs.temp ${output.dir}/g1 ${output.dir}/g2
            ""","direct_make_simple_reads_table"
        }
    }
}

direct_get_final_list = {
    doc "Get final list"
    output.dir=jaffa_output+branch
    produce(branch+".summary") {
        from(".psl", ".reads") {
            exec """
                $R --vanilla --args $input1 $input2 $transTable $knownTable $finalGapSize $exclude $output < $R_get_final_list
            ""","direct_get_final_list"
        }
    }
}

//Compile the results from multiple samples into an excel .csv table
//Make a fasta file with the candidates
direct_compile_all_results = {
    doc "Compile all results"
    var type : ""
    if (jaffa_output) {
        output.dir=jaffa_output
    }
    produce(outputName+".fasta",outputName+".csv") {
        from ("*.summary") {
            exec """
                cd ${output.dir};
                $R --vanilla --args $outputName $inputs < $R_compile_results_script ;
                function get_sequence {
                    if [ \$1 == "sample" ] ; then return ; fi ;
                    fusions_file=\$1/\$1${type}.fusions.fa ;
                    new_id=\$1---\$2---\$3 ;
                        echo ">\$new_id" >> ${outputName}.fasta ;
                    break=\$4 ;
                    sequence=`grep -A1 "^>\$3" \$fusions_file | grep -v "^>"` ;
                    start=`echo \$sequence | cut -c 1-\$((\${break}-1))` ;
                    middle=`echo \$sequence | cut -c \$break-\$((\${break}+1)) | tr '[:upper:]' '[:lower:]'` ;
                    string_length=`echo \${#sequence}` ;
                    end=`echo \$sequence | cut -c \$((\$break+2))-$string_length ` ;
                    echo ${start}${middle}${end} >> ${outputName}.fasta ;
                    `# grep \$3 \$1/\$1_genome.psl >> ${outputName}.psl ;` ;
                } ;
                rm -f ${outputName}.fasta ;
                cat ${outputName}.temp | while read line ; do get_sequence \$line ; done ;
                rm ${outputName}.temp ;
                echo "Done writing ${outputName}.fasta" ;
                echo "All Done"
            ""","direct_compile_all_results"
        }
    }
}

call_jaffa_run_check = segment { run_check }

//call_jaffa_compile_all_results = segment { compile_all_results }
//call_jaffa_compile_all_results = segment { direct_compile_all_results }

call_jaffa_direct_wrapper = segment {
    // input = sample(.1).fastq.gz
    direct_prepare_reads + // trim, fast align with bowtie2
    // input = _leftover_reads
    direct_get_unmapped +
    // output = .fasta, _discordant_pairs.bam
    // input = .fasta
    direct_align_reads_to_annotation +  // align reads to annotation using blat
    // output = .psl (blat table)
    // input = .psl
    direct_filter_transcripts +
    // output = .txt
    direct_extract_fusion_sequences +
    direct_align_transcripts_to_genome +
    direct_make_simple_reads_table +
    direct_get_final_list
}

call_jaffa_direct = segment {
    prepare_reads +

    get_unmapped + // align with bowtie2 "as unpaired", create bam with samtools, sort and index.  THEN convert fastq -> fasta, and remove duplicates

    align_reads_to_annotation +  // align reads to annotation using blat

    filter_transcripts + // R script, select reads which map to multiple genes

    extract_fusion_sequences +

    align_transcripts_to_genome +  // align fusion transcripts using blat

    make_simple_reads_table +

    get_final_list
}

// can individually override each of these JAFFA functions if need be.
// can add additional processing steps as needed.

