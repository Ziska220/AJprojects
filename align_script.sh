
 #! /usr/bin/env bash

#BSUB -J bowtie2
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -n 6

BOWTIE_IDX="/vol3/home/jhessel/ref/genomes/hg19/hg19"

FASTQ_R1="CGATGT_S1_L001_IR1_001.fastq.gz"
FASTQ_R2="CGATGT_S1_L001_IR2_001.fastq.gz"

BOWTIE_ARGS="--local"
UNALIGN_FQ="unaligned.fq.gz"
stats=align.stats.txt
bam=alignment.bam

bowtie2 $BOWTIE_ARGS -x $BOWTIE_IDX --un-conc-gz $UNALIGN_FQ \
    -1 $FASTQ_R1 -2 $FASTQ_R2 -p 6 \
    2> $stats \
    | samtools view -ShuF4 - \
    | samtools sort -o - sample.temp -m 8G \
    > $bam
samtools index $bam

CHROM_SIZE=/vol3/home/jhessel/ref/genomes/hg19/hg19.chrom.sizes
bedgraph=alignment.bg
bigwig=alignment.bw
bedtools genomecov -ibam $bam -g $CHROM_SIZE -bg > $bedgraph
bedGraphToBigWig $bedgraph $CHROM_SIZE $bigwig

