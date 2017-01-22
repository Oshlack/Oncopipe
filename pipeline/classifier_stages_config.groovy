about title: "Classifier and variant calling pipeline"

// Ensembl version 77 is equiv to GENCODE v 21
//GEN_DIR = "\$GENOMES/hg38/star/vep_compatible/77_GRCh38"
//GEN_FASTA = "\$GENOMES/hg38/vep/77_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
//GTF = "\$GENOMES/hg38/vep/77_GRCh38/Homo_sapiens.GRCh38.77.gtf"
//REF_FLAT = "/group/bioi1/rebeccae/hg38_ref_flat.txt"
//VEPDIR = "/group/bioi1/rebeccae/oncopipe/tools/ensembl-tools-release-77/scripts"
// Ensembl version 74 & 75 is equiv to GENCODE 19
// ENSEMBL fasta file has par masked on chr Y, also chr are numbered "1" instead of chr1
GEN_DIR="/group/bioi1/shared/projects/jaffa-competition/reference/star"
// lexicographically sorted
GEN_FASTA="/group/bioi1/shared/projects/jaffa-competition/reference/Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa"
GTF="/group/bioi1/shared/projects/jaffa-competition/reference/Homo_sapiens.GRCh37.75.gtf"
TRANSCRIPTS_GTF="/group/bioi1/shared/projects/jaffa-competition/reference/Homo_sapiens.GRCh37.75.no_haplotype.transcript.gtf"
//REF_FLAT = "/group/bioi1/shared/projects/jaffa-competition/reference/Homo_sapiens.GRch37.75.refflat.txt"
GATK_GEN_FASTA="\$GENOMES/hg19/1000g/human_g1k_v37.fasta" // karyotypic sorted

//genome="hg38"
//annotation="GENCODEV20"

//GEN_DIR = '$GENOMES/'+genome+'/star'
//GEN_FASTA = '$GENOMES/'+genome+'/fasta/'+genome+'.fa'
//GTF = '$GENOMES/'+genome+'/star/'+genome+'_'+annotation+'_Comp.gtf'
// don't use SAF -- or need to generate for GRCh38 ensembl 77
//SAF = '$GENOMES/'+genome+'/saf/'+genome+'_'+annotation+'_Comp.saf'

MAPPING_DIR = "mapped"

fastqc = {
    doc "Perform quality control using FASTQC"
    output.dir="fastqc"

    transform(".fastq.gz") to ("_fastqc.zip") {
        exec """
            $FASTQC -o $output.dir $inputs.gz
        ""","fastqc"
    }
    forward input
}

trim = {
    doc "Trim reads using Trimmomatic"
    output.dir="trimmed"

    filter("trim","trim") {
        exec """
            $TRIMMOMATIC/trimmomatic PE -phred33 -threads $threads
            $input1.gz $input2.gz
            $output1.gz ${output1.prefix}.unpaired.gz
            $output2.gz ${output2.prefix}.unpaired.gz
            ILLUMINACLIP:$ADAPTERS_FASTA:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:35
        ""","trim"
    }
}

star_map_1pass_PE = {
    doc "Map reads using the STAR aligner: 1st pass"
    output.dir = MAPPING_DIR+"/"+branch+"/1pass"
    def files = input1 + "," + input2

    produce ("SJ.out.tab") {
        from("*.fastq.trim.gz") {
            exec "rm -rf ${output.dir}/_STARtmp"
            exec """
                $STAR/STAR --genomeDir $GEN_DIR --readFilesIn $files --readFilesCommand zcat
                    --outSAMtype None 
                    --runThreadN $threads --outTmpDir ${output.dir}/_STARtmp/
                    --outFileNamePrefix ${output.dir}/
            ""","star"
        }
    }
}

star_gen2pass = {
    doc "Map reads using the STAR aligner: generate genome"
    output.dir = MAPPING_DIR+"/"+branch+"/genome_2pass"

    produce("Genome") {
        from ("SJ.out.tab") {
            exec "rm -rf ${output.dir}/_STARtmp"
            exec """
                $STAR/STAR --runMode genomeGenerate --genomeDir ${output.dir}
                    --genomeFastaFiles $GEN_FASTA --runThreadN $threads
                    --sjdbFileChrStartEnd ${MAPPING_DIR}/${branch}/1pass/SJ.out.tab
                    --sjdbOverhang 99 --sjdbGTFfile $GTF 
                    --outTmpDir ${output.dir}/_STARtmp/
                    --outFileNamePrefix ${output.dir}/
            ""","star"
        }
    }
}

star_map_2pass_PE = {
    doc "Map reads using the STAR aligner: 2nd pass"
    output.dir = MAPPING_DIR+"/"+branch

    // TODO: may not need Unsorted Aligned.out.bam, Aligned.toTranscriptome.out.bam
    produce ("Aligned.out.bam", "Aligned.sortedByCoord.out.bam", "Aligned.toTranscriptome.out.bam") {
        from (".fastq.trim.gz",".fastq.trim.gz") {
            exec "rm -rf ${output.dir}/_STARtmp"
            exec """
                $STAR/STAR --genomeDir ${output.dir}/genome_2pass --readFilesIn $input1,$input2
                    --readFilesCommand zcat 
                    --outSAMattrRGline ID:${branch.sample} SM:${branch.lane} LB:${sample_info[sample].library} PL:$PLATFORM PU:1
                    --outSAMtype BAM Unsorted SortedByCoordinate --runThreadN $threads
                    --quantMode TranscriptomeSAM
                    --outTmpDir ${output.dir}/_STARtmp/
                    --outFileNamePrefix ${output.dir}/
            ""","star"
        }
    }
}

// rna-seqc will only run with java 1.7 not 1.8 ...
// http://gatkforums.broadinstitute.org/dsde/discussion/8289/rna-seqc-gatk-intronicexpressionreadblock-error
rna_seqc_metrics = {
    doc "Generate RNA metrics"
    output.dir = "rna_seqc_metrics/" + branch

    produce("index.html") {
        from ("*.splitncigar.bam") {
            exec """
                $JAVA -Xmx6g -jar $RNASEQCDIR/RNA-SeQC.jar
                -s "$batch|$input|NA" -o ${output.dir}
                -r $GATK_GEN_FASTA -t $TRANSCRIPTS_GTF -ttype 2
            ""","rna_seqc"
        }
    }
}

count_reads_RNA = {
    doc "Count reads"
    output.dir = "counts"

    produce("counts.txt") {
        from ("*Aligned.sortedByCoord.out.bam") {
            exec """
                $SUBREAD/featureCounts --primary -p -t exon -g gene_name -T $threads -a $GTF -o $output $inputs
            ""","count"
//            exec """
//                $SUBREAD/featureCounts --primary -p -t exon -g GeneID -T $threads -F SAF -a $SAF -o $output $inputs
//            ""","count"
        }
    }
}

// TODO: modify to call Anthony Hawkin's ALL_Sorts
predict_class = {
    doc "Random Forest classifier"
    output.dir="classifier"

    produce("prediction.csv", "visualise.pdf") {
        from ("counts.txt") {
            exec """
                $R --vanilla --args $target_FPKM_file $target_RFmodel_file $target_topgenenames_file $inputs $output1 $output2 < $SCRIPTS/Class_Predictor.R
            ""","predict_class"
        }
    }
}

reorder_bam = {
    doc "Reorder"
    output.dir="variants"

    produce(branch+".reorder.bam") {
        from("*.dedupe.bam") {
            exec """
                $JAVA -Xmx6g -jar $PICARDDIR/picard.jar ReorderSam
                    R=$GATK_GEN_FASTA
                    I=$input.bam O=$output.bam
                    CREATE_INDEX=true
            ""","picard"
        }
    }
}

// VARIANT CALLING pipeline
// add read groups, sort
add_rg = {
    doc "Add read group"
    output.dir="variants"

    def lanes = inputs.gz.collect { (it.toString() =~ /_(L[0-9]{1,3})_/)[0][1] }.unique()
    def machines = inputs.gz.collect { (it.toString() =~ /_([A-Z0-9]*XX)_/)[0][1] }.unique()

    if (lanes.size()!=1)
        succeed report('templates/invalid_input.html') to channel: oncopipe_operator,
            subject: "Invalid input files for sample $sample: Bad lane information",
            message: """Failed to identify a unique lane number from FASTQ files: ${inputs.gz}.
                        Please check the format of the input file names""".stripIndent()
    if (machines.size()>1)
        succeed report('templates/invalid_input.html') to channel: oncopipe_operator,
            subject: "Invalid input files for sample $sample: Bad machine information",
            message: """Failed to identify a unique machine name from FASTQ files: ${inputs.gz}.
                        Please check the format of the input file names""".stripIndent()

    branch.lane = lanes[0]
    branch.machine = (machines.size() == 1) ? machines[0] : ""

    // TODO: determine library(?)
    produce(branch+".addrg.bam") {
        from ("Aligned.sortedByCoord.out.bam") {
            exec """
                $JAVA -Xmx6g -jar $PICARDDIR/picard.jar
                    AddOrReplaceReadGroups
                    I=$input.bam O=$output.bam
                    SO=coordinate
                    RGID=$sample
                    RGLB=lib1
                    RGPL=$PLATFORM
                    RGPU=${branch.lane}.${branch.machine}
                    RGSM=$sample
                    VALIDATION_STRINGENCY=LENIENT
            ""","picard"
        }
    }
}

mark_duplicates = {
    // mark duplicates, and create index
    doc "Mark duplicates and create index"
    output.dir="variants"

    produce(branch+".dedupe.bam",branch+".dedupe.metrics.txt") {
        from("*.addrg.bam") {
            exec """
                $JAVA -Xmx6g -jar $PICARDDIR/picard.jar
                    MarkDuplicates
                    I=$input.bam O=$output.bam 
                    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=$output.txt
            ""","picard"
        }
    }
}

splitncigar = {
    // split'n'trim and reassign mapping qualities
    doc "Split'n'cigar"
    output.dir="variants"

    produce(branch+".splitncigar.bam") {
        from("*.reorder.bam") {
            exec """
                $JAVA -Xmx6g -jar $GATKDIR/GenomeAnalysisTK.jar
                    -T SplitNCigarReads
                    -R $GATK_GEN_FASTA
                    -I $input.bam -o $output.bam
                    -rf ReassignOneMappingQuality
                    -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
            ""","gatk"
        }
    }
}

rnaseq_call_variants = {
    // variant calling
    doc "Call variants using GATK"
    output.dir="variants"

    var call_conf : 20.0,
        emit_conf : 20.0

    produce(branch+".gatk.vcf") {
        from("*.splitncigar.bam") {
            exec """
                $JAVA -Xmx6g -jar $GATKDIR/GenomeAnalysisTK.jar
                    -T HaplotypeCaller
                    -R $GATK_GEN_FASTA
                    -I $input.bam -o $output.vcf
                    -dontUseSoftClippedBases
                    -stand_call_conf $call_conf -stand_emit_conf $emit_conf
            ""","gatk"
        }
    }
}

filter_variants = {
    // variant filtering
    doc "Filter Variants"
    output.dir="variants"

    produce(branch.sample +".gatkfilter.vcf") {
        from("*.gatk.vcf") {
            exec """
                $JAVA -Xmx6g -jar $GATKDIR/GenomeAnalysisTK.jar
                    -T VariantFiltration
                    -R $GATK_GEN_FASTA
                    -V $input.vcf -o $output.vcf
                    -window 35 -cluster 3
                    -filterName FS -filter "FS > 30.0"
                    -filterName QD -filter "QD < 2.0"
            ""","gatk"
        }
    }
}

vep = {
    doc "Variant Effect Predictor"
    output.dir="variants"

    produce(branch.sample + ".vep.txt") {
        from("*.gatkfilter.vcf") {
            // requires module perl/5.20.3 on meerkat
            exec """
                perl $VEPDIR/variant_effect_predictor/variant_effect_predictor.pl
                    -i $input.vcf -o $output.txt --force_overwrite
                    --cache --dir $VEPCACHE
                    -fork $threads --symbol
                    --check_alleles --filter_common -sift b -polyphen b
            ""","vep"
        }
    }
}

vepfilter = {
    doc "Filter Variants"
    output.dir="variants"

    produce(branch.sample + ".vepfilter.txt") {
        from("*.vep.txt") {
            exec """
                perl $VEPDIR/variant_effect_predictor/filter_vep.pl
                    -i $input.txt -o $output.txt --force_overwrite
                    --filter "SYMBOL in $target_gene_file"
            ""","vep_filter"
        }
    }
}

vepvcf = {
    doc "VEP VCF"
    output.dir="variants"

    produce(branch.sample + ".vep.vcf") {
        from("*.gatkfilter.vcf") {
            exec """
                perl $VEPDIR/variant_effect_predictor/variant_effect_predictor.pl
                    -i $input.vcf -o $output.vcf --force_overwrite
                    --cache --dir $VEPCACHE
                    --vcf --allele_number -fork $threads --symbol
                    --check_alleles --filter_common -sift b -polyphen b
            ""","vep"
        }
    }
}

vepfiltervcf = {
    doc "Filter VCF"
    output.dir="variants"

    produce(branch.sample +".vepfilter.vcf") {
        from("*.vep.vcf") {
            exec """
                perl $VEPDIR/variant_effect_predictor/filter-vep.pl
                    -i $input.vcf -o $output.vcf --force_overwrite
                    --filter "SYMBOL in $target_gene_file"
            ""","vep_filter"
        }
    }
}

// Picard CollectRnaSeqMetrics metrics
collectRNAMetrics = {
    doc "Collect RNA metrics"
    output.dir="variants"

    produce("output.RNA_Metrics") {
        from("*.addrg.bam") {
            exec """
                $JAVA -jar $PICARDDIR/picard.jar CollectRnaSeqMetrics I=$input.bam 0=$output.RNA_Metrics REF_FLAT=$REF_FLAT STRAND=NONE
            ""","picard"
        }
    }
}

classifier_pipeline = segment {
    trim + 

    // align to genome
    star_map_1pass_PE +
    star_gen2pass +
    star_map_2pass_PE +

    // variant calling pipeline (Picard, GATK)
    add_rg +
    mark_duplicates +
    reorder_bam +
    splitncigar +
    rna_seqc_metrics +
    rnaseq_call_variants +
    filter_variants //+
//    vepvcf +
//    vepfiltervcf +
//    vep
}
