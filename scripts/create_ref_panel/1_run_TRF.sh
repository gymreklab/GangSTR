#!/bin/bash


REF=PATH/TO/CHROM/SEPARATED/REF/ # Folder containing chrom separated reference genome (chr1.fa chr2.fa chr3.fa ...)
OUT=PATH/TO/OUT.bed # Output file
TRF=trf # Path to TRF executable (if trf is not added to path)
THREADS=4 # Number of threads
MAX_CHR_NUM=22 # Maximum chr number (other than sex chroms) human:22, mouse:19
NON_NUM_CHRS="X,Y"
TMP=PATH/TO/TMP/ # Path to folder to store temp files

scripts/chr_make_reference.sh \
    $REF \
    $OUT \
    trf \
    4 \
    $MAX_CHR_NUM \
    $NON_NUM_CHR \
    $TMP
