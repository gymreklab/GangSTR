#!/bin/bash

./src/GangSTR \
    --bam test.bam \
    --ref /storage/resources/dbase/human/hs37d5/hs37d5.fa \
    --regions test.bed \
    --out test
