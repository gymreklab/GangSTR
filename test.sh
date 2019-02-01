#!/bin/bash

./src/GangSTR \
    --bam /storage/nmmsv/plat_genome/NA12878_S1.bam \
    --ref /storage/resources/dbase/human/hg19/hg19.fa \
    --regions /storage/resources/dbase/human/hg19/pathogenic_loci_hg19.bed \
    --out test --verbose
