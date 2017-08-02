#!/bin/bash
make -j -C ../

locus=HTT

bam_dir=/storage/resources/datasets/repeat-expansions/bams

bed_dir=../tests/$locus.bed

for bam_file in $bam_dir/*.bam;
# for bam_file in $bam_dir/SRR4243223_47_20.rmdup.bam;
do
	echo '>>Running for: '$bam_file
	../src/GangSTR \
	    --bam $bam_file \
	    --ref /storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta \
	    --regions $bed_dir \
	    --out test
done

