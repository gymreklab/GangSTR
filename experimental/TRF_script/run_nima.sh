#!/bin/bash
#hg19: /storage/resources/dbase/human/hg19/chromFa/ \
# New hg19: /storage/nmmsv/ref_genome/hg19/chromFa/ \
#hg38: /storage/nmmsv/ref_genome/chromFa/
# mm10: /storage/nmmsv/ref_genome/mouse/mm10
./chr_make_reference.sh \
    /storage/nmmsv/ref_genome/mouse/mm10/ \
    /storage/nmmsv/reference_TRF/raw_mm10_20bp_XY_mis3_ins5.bed \
    trf \
    4 \
    19 #maximum chr number (other than sex chroms) human:22, mouse:19
