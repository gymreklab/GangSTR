#!/bin/bash

BED=/storage/nmmsv/reference_TRF/hg19/hg19_ver13.bed
DISEASE=/storage/nmmsv/reference_TRF/disease/hg19_disease_thresh.bed
READLEN=150
OUT=/storage/nmmsv/analysis/GangSTR-analyses/genome_wide/str_info/hg19_ver13_strinfo_pop_pad3.bed
POPULATION=/storage/resources/datasets/ENA-illuminia-150/ENA-illumina-150.strinfo.hg19.bed

python info_gen.py \
    --bed $BED \
    --out $OUT \
    --readlen $READLEN \
    --disease $DISEASE \
    --population $POPULATION \
    --population-pad 3 \
    


