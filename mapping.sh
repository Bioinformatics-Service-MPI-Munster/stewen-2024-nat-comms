OUTPUT_FOLDER=/home/kkruse/project-results/stewen/huvec
FASTQ1_FILE=/home/research/data/scRNA-Seq_raw/20200626_AAJ_HUVEC-HUAEC/Sample-3_S3_R1_001.fastq.gz
FASTQ2_FILE=/home/research/data/scRNA-Seq_raw/20200626_AAJ_HUVEC-HUAEC/Sample-3_S3_R2_001.fastq.gz
WHITELIST_FILE=/home/kkruse/data/rhapsody/rhapsody_whitelist.txt
INDEX=/home/kkruse/data/index/star/hg-grch38.p13

mkdir "${OUTPUT_FOLDER}/trimmed_fastq"
trim_galore \
  --length 66 \
  --paired \
  --gzip \
  -q 20 \
  -o "${OUTPUT_FOLDER}/trimmed_fastq" \
  $FASTQ1_FILE \
  $FASTQ2_FILE

mkdir "${OUTPUT_FOLDER}/fastq"
rhapsody-extract-barcode \
  -w $WHITELIST_FILE \
  -b ${OUTPUT_FOLDER}/fastq/barcodes.txt \
  ${OUTPUT_FOLDER}/trimmed_fastq/Sample-3_S3_R1_001_val_1.fq.gz \
  ${OUTPUT_FOLDER}/trimmed_fastq/Sample-3_S3_R2_001_val_2.fq.gz \
  ${OUTPUT_FOLDER}/fastq/stewen_huvec_R1_barcodes.fastq.gz \
  ${OUTPUT_FOLDER}/fastq/stewen_huvec_R2_reads.fastq.gz

mkdir "${OUTPUT_FOLDER}/fastq/dmx/"
rhapsody-demultiplex \
    ${OUTPUT_FOLDER}/fastq/stewen_huvec_R1_barcodes.fastq.gz \
    ${OUTPUT_FOLDER}/fastq/stewen_huvec_R2_reads.fastq.gz \
    hg \
    ${OUTPUT_FOLDER}/fastq/dmx/ \
    -p stewen_huvec \
    -n ${OUTPUT_FOLDER}/fastq/dmx/stewen_huvec_noise.png \
    -c ${OUTPUT_FOLDER}/fastq/dmx/stewen_huvec_sample_tag_counts.png


STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ${INDEX} \
  --genomeFastaFiles ~/data/genomes/GRCh38.p13.genome.fa \
  --sjdbGTFfile ~/data/genome-annotations/gencode/hg/gencode.v35.annotation.gtf


run_starsolo () {
    STAR \
        --runThreadN 16 \
        --genomeDir ${INDEX} \
        --soloType CB_UMI_Simple \
        --soloCellFilter None \
        --outSAMtype BAM SortedByCoordinate \
        --soloCBwhitelist ${WHITELIST_FILE} \
        --soloFeatures Gene  \
        --soloCBstart 1 \
        --soloCBlen 27 \
        --soloUMIlen 8 \
        --soloUMIstart 28 \
        --readFilesCommand zcat \
        --readFilesIn ${READS_FILE} ${BARCODES_FILE} \
        --outFileNamePrefix ${OUTPUT_PREFIX}
}

INDEX=/home/kkruse/data/index/star/hg-grch38.p13

READS_FILE=${OUTPUT_FOLDER}/fastq/dmx/stewen_huvec_5_R2.fastq.gz
BARCODES_FILE=${OUTPUT_FOLDER}/fastq/dmx/stewen_huvec_5_R1.fastq.gz
OUTPUT_PREFIX=${OUTPUT_FOLDER}/star/stewen_huvec_si_ctrl/
mkdir ${OUTPUT_PREFIX}
run_starsolo

READS_FILE=${OUTPUT_FOLDER}/fastq/dmx/stewen_huvec_6_R2.fastq.gz
BARCODES_FILE=${OUTPUT_FOLDER}/fastq/dmx/stewen_huvec_6_R1.fastq.gz
OUTPUT_PREFIX=${OUTPUT_FOLDER}/star/stewen_huvec_si_ephb4/
mkdir ${OUTPUT_PREFIX}
run_starsolo

READS_FILE=${OUTPUT_FOLDER}/fastq/dmx/stewen_huvec_7_R2.fastq.gz
BARCODES_FILE=${OUTPUT_FOLDER}/fastq/dmx/stewen_huvec_7_R1.fastq.gz
OUTPUT_PREFIX=${OUTPUT_FOLDER}/star/stewen_huvec_si_efnb2/
mkdir ${OUTPUT_PREFIX}
run_starsolo

