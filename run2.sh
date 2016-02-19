#! /usr/bin/env bash

datasets='/Users/test/GitHub/ProblemSet1/problem-set-1/problem-set-2/data-sets'
#Use BEDtools intersect to identify the size of the largest overlap
#between CTCF and H3K4me3 locations.

tfbs="$datasets/bed/encode.tfbs.chr22.bed.gz"
h3k4me3="$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz"

answer_1=$(bedtools intersect -a $tfbs -b $h3k4me3 -wo \
    | awk '($4 == "CTCF") {print$3-$2}' \
    | sort -k1nr \
    | head -n1)

echo "answer_1: $answer_1" > answers.yml

#Use BEDtools to calculate the GC content of nucleotides 19,000,000 to
#19,000,500 on chr22 of hg19 genome build. Report the GC content as a
#fraction (e.g., 0.50)

chr22fafile="$datasets/fasta/hg19.chr22.fa"

#answer_2=echo -e "chr22\t19000000\t190005000" > posnuc.bed
echo -e "chr22\t19000000\t19000500" > posnuc.bed

a2=$(bedtools nuc -fi $chr22fafile -bed posnuc.bed | cut -f5 | tail -n1) 

echo "answer_2: $a2" >> answers.yml


#Use BEDtools to identify the length of the CTCF ChIP-seq peak (i.e.,
#interval) that has the largest mean signal in ctcf.hela.chr22.bg.gz.

peakschr22="$datasets/bed/peaks.chr22.bed.gz"
ctcfhelachr22="$datasets/bedtools/ctcf.hela.chr22.bg.gz"

answer_3=$(bedtools map -a $peakschr22 -b $ctcfhelachr22 -c 4 -o mean \
    | awk '($4 == "CTCF")' \
    | sort -k5nr \
    | awk '{print $3-$2}' \
    | head -n1)

echo "answer_3: $answer_3" >> answers.yml

#Use BEDtools to identify the gene promoter (defined as 1000 bp upstream
#of a TSS) with the highest median signal in ctcf.hela.chr22.bg.gz. Report
#the gene name (e.g., 'ABC123')

tss="$datasets/bed/tss.hg19.chr22.bed.gz"
hg19="$datasets/genome/hg19.genome"
ctcf="$datasets/bedtools/ctcf.hela.chr22.bg.gz"

answer_4=$(bedtools flank -b 1000 -i $tss -g $hg19 | bedtools sort -i - \
    | bedtools map -a - -b $ctcf -c 4 -o median \
    | awk '$7 != "."' \
    | sort -k7nr \
    | head -n1 \
    | cut -f4)

echo "answer_4: $answer_4" >> answers.yml

#Use BEDtools to identify the longest interval on chr22 that is not
#covered by genes.hg19.bed.gz. Report the interval like chr1:100-500.

genes="$datasets/bed/genes.hg19.bed.gz"
hg19genome="$datasets/genome/hg19.genome"

answer_5=$(bedtools complement -i $genes -g $hg19genome \
    | bedtools sort -i - \
    | awk 'BEGIN{OFS ="\t"} {print $0, ($3-$2)}' \
    | sort -k4nr \
    | awk '{print $1, ":", $2, "-", $3}' \
    | head -n1)

echo "answer_5: $answer_5" >> answers.yml


