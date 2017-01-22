////////////////////////////////////////////////////////////////////////
//
// This is the main configuration file for oncopipe.
//
// In here are configured the locations for all the tools that the
// pipeline uses. You should copy this file to 'config.groovy'
// and then read through instructions below to edit it for your 
// local setup. 
//
// NOTE: please use C-style single line comments (//), and not BASH 
// style comments in this file (#)
//
////////////////////////////////////////////////////////////////////////
/////////////////////////// BASIC PARAMETERS ///////////////////////////
//
// The base of everything - set this to the absolute path of the 
// root of the pipeline distribution (most likely, parent folder of
// the folder this file is in)
BASE="/group/bioi1/rebeccae/oncopipe"

// Set a good location for storing large temp files here (probably not /tmp)
TMPDIR="$BASE/tmpdata"

// Enter email here to get notified by email about failures
EMAILS=""

// If you are using the default reference data, download it now in 
// the hg19 folder (see hg19/README)
//
// NEXT: run ./pipeline/scripts/install.sh from 
// the root of the distribution

//////////////////// REFERENCE FILES ////////////////////////////////////
//
// You do not need to edit below here if you are using the default 
// HG19 reference files. These can be downloaded by the installer script.
// However if you want to avoid downloading those or use your own reference
// files then you should enter correct paths below.
//
////////////////////////////////////////////////////////////////////////

// Set location of your reference files here (see hg19/README for what is required)
REFBASE="$BASE/hg19" 

// Set to the reference FASTA file, which must be indexed with bwa, and samtools faidx
// To download this file, and the files below, visit the GATK resource bundle
// at ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19
REF="$REFBASE/ucsc.hg19.fasta"

// Set to a VCF file containing DBSNP entries (or leave it if you are downloading the default)
DBSNP="$REFBASE/dbsnp_138.hg19.vcf"

// Set to a VCF file containing known indels here
GOLD_STANDARD_INDELS="$REFBASE/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"

// For self tests and other default features to work, you should
// set a "default" exome target here. Note that you can specify a different
// exome capture region for any individual analysis
EXOME_TARGET=""

//GEN_DIR = '$GENOMES/hg38/star'
//GEN_FASTA = '$GENOMES/hg38/fasta/hg38.fa'
//GTF = '$GENOMES/hg38/star/hg38_GENCODEV20_comp.gtf'
//SAF = '$GENOMES/hg38/saf/hg38_GENCODEV20_Comp.saf'


//////////////////// OPTIONAL PARAMTERS ////////////////////////////////
//
// You probably do NOT need to set anything under here! 
// 
// However you may wish to customise the locations of some tools such
// as java, python and perl installations.
// 
///////////////////////////////////////////////////////////////////////

// The variant database requires every sample id to map to a different
// individual. If multiple sample ids can map to the same individual
// (for example, you repeat sequencing on a sample, etc.), then
// you can mask out the part of the sample id that is not unique
// with a regular expression here. For example, Melbourne Genomics
// uses a 9 digit sample (study) id, but only the first 7 digits are
// unique for an individual. We can use a mask of ".{7}" to indicate
// that only the first 8 digits of the study id should be used in
// the variant database.
SAMPLE_ID_MASK=".*"

// Filter out variants observed more than 10 times (ie: 11 times or more)
// in samples from a different cohort / disease target
OUT_OF_COHORT_VARIANT_COUNT_FILTER=10

// This is only used for setting read group information in the
// BAM files
PLATFORM="illumina"

// The coverage level below which a sample should be failed
MEDIAN_COVERAGE_THRESHOLD=50

// Base location of all the tools that we use
TOOLS="$BASE/tools"

// Various support scripts that the pipeline uses
SCRIPTS="$BASE/pipeline/scripts"

// Location of Picard tools here
PICARD_HOME="$BASE/tools/picard/picard-tools-1.65"

// Set location of Annovar distribution
// Due to license restrictions on Annovar, you must download
// it yourself and place it in this location
// Note also that many databases also need to be downloaded,
// using annovar's downdb function. See scripts/download_annovar_db.sh for
// a helper script.
ANNOVAR="$TOOLS/annovar/2015-03-22"
ANNOVAR_DB="$TOOLS/annovar/humandb"

JAFFA="$TOOLS/JAFFA"

// Due to GATK license restrictions, you must download
// GATK yourself and place it in this location
// (or point this to your installation)
// GATK="$TOOLS/gatk/2.3.9"
GATK="$TOOLS/gatk/2.3.9"
GATK_LEGACY=true

// Utilities for making Excel files 
EXCEL="$TOOLS/excel/1.0"

// Location of Bedtools distribution
// dev.meerkat
BEDTOOLS="$TOOLS/bedtools/2.25.0"
// mbioapp1
//BEDTOOLS="$TOOLS/bedtools/2.18.2"

// Location of Samtools
SAMTOOLS="$TOOLS/samtools"

// FastQC tool
FASTQC="$TOOLS/fastqc"

// Utilties for processing NGS data in Java/Groovy
GROOVY_NGS="$TOOLS/groovy-ngs-utils/1.0.2"

// Set location of Variant Effect Predictor here
// and store it in the local directory called 'vep_cache'
// (you can create a symlink to an existing directory with 
// that name if desired).
// See tools/vep/README for more information
//VEP_VERSION="74"
//VEP="$TOOLS/vep/$VEP_VERSION"
VEPCACHE="/group/bioi1/rebeccae/cpipe/tools/vep/vep_cache/"

IGVTOOLS="$TOOLS/IGVTools/2.3.6"

// IGV location
IGV="$TOOLS/tools/igv/2.3.15"

// Location and version of BWA
BWA="$TOOLS/bwa/0.7.5a/bwa"
BWA_THREADS="5"

CONDEL="$TOOLS/condel/24-1-2013"

// Database of unique variants, updated for each sample
VARIANT_DB="$BASE/variants.db"
ID_FILE="$BASE/pipeline_id"

// By default variant counts are annotated from the same database as the
// one that they were added to in the first place
// However you can modify them to be separate if you wish
UPDATE_VARIANT_DB=VARIANT_DB
ANNOTATION_VARIANT_DB=VARIANT_DB

// Location of groovy installation
GROOVY_HOME="$TOOLS/groovy/2.3.4"

// GROOVY binary
GROOVY="$GROOVY_HOME/bin/groovy"

// Whether to fail analysis if FASTQC produces warnings
CHECK_FASTQC_FAILURES=false

// If adapter sequence for exome capture is known, put here
ADAPTERS_FASTA=false

// By default synonymous variants are excluded from all outputs
// Set to true to include them
EXCLUDE_VARIANT_TYPES="synonymous SNV"

// Use default Java installed in PATH
JAVA="java"

// SNPEFF location
// Note: snpeff is not needed by default pipeline
SNPEFF="$TOOLS/snpeff/3.1"

// Genome needed for expanded splice regions
HG19_CHROM_INFO="$REFBASE/hg19.genome"

// Trimmomatic location
TRIMMOMATIC="$TOOLS/trimmomatic"
ADAPTERS_FASTA="/usr/local/installed/trimmomatic/0.35/adapters/TruSeq2-PE.fa"

STAR="$TOOLS/star"
SUBREAD="$TOOLS/subread"

// Bamsurgeon location
BAMSURGEON="$TOOLS/bamsurgeon/20150331"


PYTHON="python"

splice_region_window=2

// interval padding to pass to the variant caller
INTERVAL_PADDING_CALL=25

// interval padding for SNVs in filter_variants
INTERVAL_PADDING_SNV=10

// interval padding for indels in filter_variants
INTERVAL_PADDING_INDEL=25

// do not filter synonymous with this range
ALLOW_SYNONYMOUS_INTRON=0
ALLOW_SYNONYMOUS_EXON=0
