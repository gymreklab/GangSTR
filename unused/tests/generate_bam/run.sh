#!/bin/bash

set -e

bam=../test.enclosing.single.bam
ref=/storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta

/storage/nmmsv/GangSTR/src/GangSTR \
    --stutterup 0.001 --stutterdown 0.001 --stutterprob 0.99 \
    --frrweight 0.5 --enclweight 1.0 --spanweight 1.0 --flankweight 1.0 \
    --bam $bam \
    --ref $ref \
    --regions test_locus.bed \
    --ploidy 2 \
    --minmatch 5\
    --insertmax 2000\
    --insertmean 500\
    --insertsdev 50\
    --output-bootstraps --output-readinfo --out test --verbose

bgzip -f test.vcf
tabix -f -p vcf test.vcf.gz
bcftools query -f '%CHROM\t%POS\t%END\t[%GB\t%CI\t%RC]\n' test.vcf.gz
