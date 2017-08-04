#!/bin/bash

locus=ATXN3
exp_name=cpp_ATXN3_8_cov60_dist200_sd70_rl76_DIP_const30
const=30
nc_list='1 3 5 7 10 12 15 20 25 30 40 50 60 70 80 90 100 120 150'
echo $exp_name
for nc in $nc_list;
do
	echo '>>Running for: ('$nc,$const')'
	./src/GangSTR \
	--bam ../expansion-experiments/$exp_name/aligned_read/nc_$nc.sorted.bam \
	--ref /storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta \
	--regions tests/$locus.bed \
	--out test
done

