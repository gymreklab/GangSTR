#!/bin/bash

./src/GangSTR \
    --bam /storage/resources/datasets/gtex/bams/SRR2157423.bam \
    --ref /storage/resources/dbase/human/hs37d5/hs37d5.fa \
    --regions test.bed \
    --out test
