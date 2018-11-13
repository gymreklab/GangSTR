#!/bin/bash
#hg19: /storage/resources/dbase/human/hg19/chromFa/ \
#hg38: /storage/nmmsv/ref_genome/chromFa/
./chr_make_reference.sh \
    /storage/resources/dbase/human/hg19/chromFa/ \
    /storage/nmmsv/reference_TRF/raw_ref_hg19_20bp_XY.bed \
    trf \
    4
