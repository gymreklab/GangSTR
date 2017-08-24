#!/bin/bash

locus=HTT

bam_dir=/storage/resources/datasets/repeat-expansions/bams
bam_id=SRR4243174_43_22.rmdup.bam
read_id=SRR4243174.9149378


# samtools view $bam_dir/$bam_id 4:3063003-3089903 | grep $read_id | less -S
# samtools view $bam_dir/$bam_id 18:53233385-53273385 | grep $read_id | less -S
samtools view $bam_dir/$bam_id 4:189142085-189222085 | grep $read_id | less -S
