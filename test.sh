#!/bin/bash

./src/GangSTR \
    --bam /storage/resources/datasets/PlatinumGenomesDbGaP/sra/SRR4435260.bam,/storage/resources/datasets/PlatinumGenomesDbGaP/sra/SRR4435261.bam \
    --ref /storage/resources/dbase/human/hg19/hg19.fa \
    --regions test.bed \
    --out test --verbose
