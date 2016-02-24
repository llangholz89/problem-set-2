#! usr/bin/env bash

#Problem Set 2 - Bedtools

#Files:
fasta_chr1='../data-sets/bedtools/chr1.hg19.fa.gz'
histone_bed='../data-sets/bed/encode.h3k4me3.hela.chr22.bed.gz'
hg19genes='../data-sets/bed/genes.hg19.bed.gz'
tfChIP='../data-sets/bed/encode.tfbs.chr22.bed.gz'
ctcf='../data-sets/bedtools/ctcf.hela.chr22.bg.gz'
genome='../data-sets/bedtools/hg19.genome'
tss='../data-sets/bed/tss.hg19.chr22.bed.gz'
fasta_chr22='../data-sets/fasta/hg19.chr22.fa'

#----------------------------------------------
#1) Identify size of the largest overlap between CTCF and H3K4me locations

#find ctcf peaks and put in bed file:
gzcat $tfChIP | awk '$4 == "CTCF"' > ctcf-peaks.bed
#intersect ctcf peaks with H3K4me peaks
answer1=$(bedtools intersect -a ctcf-peaks.bed -b $histone_bed -wo \
    | awk '{print $NF}' \
    | sort -nr \
    | head -n1) 
#note: "-wo" puts a column at the end with the number of bases overlapping
echo "answer-1: $answer1" > answers.yml

#---------------------------------------------
#2) Calculate the GC content of nucleotides 19,000,000 to 19,000,500 on
#chr22 of hg19 genome and report as a fraction

#create a bed file:
echo -e "chr22\t19000000\t19000500" > q2-new.bed
answer2=$(bedtools nuc -fi $fasta_chr22 -bed q2-new.bed \
    | grep -v '^#' \
    | awk '{print $5}')

echo "answer-2: $answer2" >> answers.yml

#----------------------------------------------
#3) length of the CTCF ChIP-seq peak (interval) with largest mean signal
answer3=$(bedtools map -a ctcf-peaks.bed -b $ctcf -c 4 -o mean \
    |sort -k5rn \
    |head -n1 \
    |awk '{print $3-$2}')

echo "answer-3: $answer3" >> answers.yml

#----------------------------------------------
#4) Identify gene promoter (1000bp upstream of TSS) with highest median
# signal and report gene name

#create promoter intervals:
#bedtools flank -i $tss -g $genome -b 1000 | sort -k2n > ctcf-promoters.bed

#map signal onto promoter intervals and obtain median:
answer4=$(bedtools map -a ctcf-promoters.bed -b $ctcf -c 4 -o median \
    | sort -k7rn \
    | head -n1 \
    | cut -f4)

echo "answer-4: $answer4" >> answers.yml

#----------------------------------------
#5) Identify longet interval on chr22 not covered by genes.hg19.bed.gz
answer5=$(bedtools complement -i $hg19genes -g $genome \
    | awk '($1=="chr22") {print $1, $2, $3, $3-$2}' \
    | sort -k4nr \
    | head -n1 \
    | awk '{print $1":"$2"-"$3}')

echo "answer-5: $answer5" >> answers.yml

#----------------------------------------
#6) Extra Credit
#using the bedtools "links" function:
bedtools links -i ctcf-peaks.bed > ctcf-peaks.html
# Creates an html file with links to the UCSC Genome Browser for all of
# our peak intervals in the ctcf-peaks bed file and allows for manual
# inspection of these intervals 









