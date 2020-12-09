#!/bin/bash

######################################################################################################
# File name: RNA_seq_pipeline.sh
# Author: Huairen Zhang
# Institution: Institute of genetics and developmental biology, Chinese Academic Sciences
# Email: samhuairen@gmail.com
# Created Time: Thu 6 Oct. 2020 17:45 PM
######################################################################################################

if [[ "$1" == "-h" || "$1" == "-help" || "$1" == "--help" ]];
then
    echo "##################################################################################"
    echo "Performs reads alignment,  select unique maping reads for downstream analysis,
    describes mapping staticstics,converts bam files and bw files"
    echo ""
    echo "Require hisat2, samtools, read_distribution.py,bam2bigwig.py, HTSeq-count.py "
    echo ""
    echo "Usage: RNA_seq_pipeline.sh <*.fastq1 file> <*.fastq2 file> <number of base without mapping from 5UTR>
     <number of base without mapping from 3UTR> <number of threads> <*.fasta reference file>
      <PATH to hisat2 mapping index> <*.gene annotation bed file> <gff3 file> <HT-Seq count.py path> <file prefix>"
    echo "#########################################################################################"
    exit 0
fi

filein1="$1" # Paired-end forward *.fastq1 file
filein2="$2" # Paired-end reverse *.fastq2 file
num_5utr_cut="$3" # Number of base without mapping from 5'URT
num_3utr_cut="$4" # Number of base without mapping from 3'URT
n_treads="$5" # Number of threads used for hisat2 mapping and samtools
genome_fasta="$6" # Reference genome fasta file
index_path="$7"  # Hisat2 align index folder
bed_file_path="$8" # Bed file from whole genome annotation protein
gff3="$9" #gff3 annotation_file
HTSEQ_COUNT="${10}" #HT_Seq count.py path 
prefix="${11}" # file prefix
##################################################################################################
echo "#################################################################################################################"
echo ""
echo "Hisat2 mapping starts..."
hisat2  -5 $num_5utr_cut -3 $num_3utr_cut -p $n_treads  -x $index_path -1 $filein1 -2  $filein2 -S $prefix.sam  2>$prefix.hisat2.log
echo "hisat2 mapping complete!"
echo ""

#################################################################################################
echo "extract unique reads mapped on the genome for the downstream analysis"
grep "NH:i:1" $prefix.sam  > $prefix.unique.sam
samtools view -H   $prefix.sam > $prefix.header
cat $prefix.header $prefix.unique.sam > $prefix.unique.aligned.sam
echo "extration unique mapped reads complete!"
echo ""


###################################################################################################
echo "samtools convert unique sam file to bam file"
samtools view -@ $n_treads  -S -b $prefix.unique.aligned.sam > $prefix.unique.aligned.bam
samtools sort -@ $n_treads   $prefix.unique.aligned.bam -o $prefix.unique.aligned_sorted.bam
samtools index -@ $n_treads  $prefix.unique.aligned_sorted.bam
echo "sam to bam conversion complete!"
echo ""

###################################################################################################
echo "HTSeq_count to count reads on feature"
python $HTSEQ_COUNT $prefix.unique.aligned_sorted.bam $gff3 -f bam -r pos -s no -t gene -i ID 1>${prefix}.geneCounts 2>${prefix}.HTseq.log
echo "done"
echo ""


####################################################################################################
echo "using reads_distribution.py computer reads distribution across UTR,CDS, TSS up and down region"
read_distribution.py -i $prefix.unique.aligned_sorted.bam -r $bed_file_path > $prefix.reads.distribution.txt
echo "Read_distribution calculation done!"
echo ""


####################################################################################################
echo "convert bam files to bw files for viewing with IGV"
bam2bigwig.py -i $prefix.unique.aligned_sorted.bam  -g $genome_fasta  -s both -o $prefix.bw
echo "bam2bigwin complete!"
echo ""


#####################################################################################################
echo "delete middle sam and bam files"
rm $prefix.sam
rm $prefix.unique.sam
rm $prefix.header
rm $prefix.unique.aligned.sam
rm $prefix.unique.aligned.bam
echo "done!"
echo ""
