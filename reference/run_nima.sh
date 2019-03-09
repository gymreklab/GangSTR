#!/bin/bash
#hg19: /storage/resources/dbase/human/hg19/chromFa/ \
# New hg19: /storage/nmmsv/ref_genome/hg19/chromFa/ \
#hg38: /storage/nmmsv/ref_genome/chromFa/
./chr_make_reference.sh \
    /storage/nmmsv/ref_genome/hg19/chromFa/ \
    /storage/nmmsv/reference_TRF/raw_ref_hg19_20bp_XY_mis3_ins5.bed \
    trf \
    4
