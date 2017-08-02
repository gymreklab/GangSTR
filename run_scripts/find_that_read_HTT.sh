#!/bin/bash

locus=HTT

bam_dir=/storage/resources/datasets/repeat-expansions/bams
bam_id=SRR4243189_44_23.rmdup.bam
read_id=SRR4243189.38822722


samtools view $bam_dir/$bam_id 4:3074003-3079903 | grep $read_id | less -S