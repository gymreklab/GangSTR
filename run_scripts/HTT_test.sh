#!/bin/bash
repo_dir=/storage/nmmsv/GangSTR/
make -j -C $repo_dir

locus=HTT

bam_dir=/storage/resources/datasets/repeat-expansions/bams


bed_dir=$repo_dir/tests/$locus.bed

# for bam_file in $bam_dir/*.bam;
for bam_file in $bam_dir/SRR4243174_43_22.rmdup.bam;
do
	echo '>>Running for: '$bam_file
	$repo_dir/src/GangSTR \
	    --bam $bam_file \
	    --ref /storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta \
	    --regions $bed_dir \
	    --out test \
	    --frrweight 0.12 \
	    --enclweight 0.3 \
	    --spanweight 1.0 \
	    --flankweight 0.5 \
	    --ploidy 2 \

done

