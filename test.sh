#!/bin/bash

./src/GangSTR \
    --bam /storage/resources/datasets/PlatinumGenomesDbGaP/sra/SRR4435259.bam \
    --ref /storage/resources/dbase/human/hg19/hg19.fa \
    --regions /storage/resources/dbase/human/hg19/pathogenic_loci_hg19.bed \
    --out test --verbose --include-ggl --chrom chr19
