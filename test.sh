#!/bin/bash

./src/GangSTR \
    --bam tests/test.sorted.bam \
    --ref /storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta \
    --regions test.bed \
    --out test
