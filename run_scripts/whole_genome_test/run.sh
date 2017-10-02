#!/bin/bash

set -e

NA12878=/storage/resources/datasets/GM12878/platinum_genomes_200x/bam/NA12878_platinum_genomes.bam
ref=/storage/resources/dbase/human/hg19/Homo_sapiens_assembly19.fasta

ref_chr=/storage/resources/dbase/human/hg19/hg19.fa
# SBMA (ChrX) Expanded:
NA23709=/storage/resources/datasets/IlluminaRepeatExpansions/forExpValidation/_EGAR00001587365_repeat_expansions_NA23709.bam

# ATXN3 SCA3 (Chr14) Expanded:
NA06151=/storage/resources/datasets/IlluminaRepeatExpansions/forExpValidation/_EGAR00001587277_repeat_expansions_NA06151.bam

# ATXN1 SCA1 (Chr6)
NA06926=/storage/resources/datasets/IlluminaRepeatExpansions/forExpValidation/_EGAR00001587302_repeat_expansions_NA06926.bam

# HTT (Chr4) Expanded
NA13509=/storage/resources/datasets/IlluminaRepeatExpansions/forExpValidation/_EGAR00001587344_repeat_expansions_NA13509.bam


/storage/nmmsv/GangSTR/src/GangSTR \
    --stutterup 0.001 --stutterdown 0.001 --stutterprob 0.99 \
    --frrweight 0.5 --enclweight 1.0 --spanweight 1.0 --flankweight 1.0 \
    --bam $NA06151 \
    --ref $ref_chr \
    --regions /storage/nmmsv/GangSTR/reference/test_locus.bed \
    --ploidy 2 \
    --minmatch 5\
    --output-bootstraps --output-readinfo --out NA12878_validation --verbose

bgzip -f NA12878_validation.vcf
tabix -f -p vcf NA12878_validation.vcf.gz
bcftools query -f '%CHROM\t%POS\t%END\t[%GB\t%CI\t%RC]\n' NA12878_validation.vcf.gz
