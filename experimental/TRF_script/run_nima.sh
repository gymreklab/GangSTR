#!/bin/bash
#hg19: /storage/resources/dbase/human/hg19/chromFa/ \
# New hg19: /storage/nmmsv/ref_genome/hg19/chromFa/ \
#hg38: /storage/nmmsv/ref_genome/chromFa/
# mm10: /storage/nmmsv/ref_genome/mouse/mm10

#maximum chr number (other than sex chroms) human:22, mouse:19

./chr_make_reference.sh \
    /storage/nmmsv/ref_genome/mouse/mm10/ \
    /storage/nmmsv/reference_TRF/raw_mm10_20bp_XY.bed \
    trf \
    4 \
    19 \
    /storage/nmmsv/tmp/
