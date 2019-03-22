#!/bin/bash

BED=/storage/nmmsv/reference_TRF/hg19/hg19_ver13.bed
OVERWRITE=/storage/nmmsv/reference_TRF/disease/hg19_disease_thresh.bed
READLEN=150
OUT=/storage/nmmsv/analysis/GangSTR-analyses/genome_wide/str_info/hg19_ver13_strinfo_dis.bed


python str-info-generator.py \
    --bed $BED \
    --out $OUT \
    --readlen $READLEN \
    --overwrite $OVERWRITE
    


