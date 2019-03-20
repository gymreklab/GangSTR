#!/bin/bash


VER=15
REF="hg38"
if [ "$REF" = "hg19" ]; then
    TRF_IN=/storage/nmmsv/reference_TRF/raw/raw_ref_hg19_20bp_XY.bed
    OUT=/storage/nmmsv/reference_TRF/hg19_ver${VER}_tmp.bed
elif [ "$REF" = "hg38" ]; then
    TRF_IN=/storage/nmmsv/reference_TRF/raw/raw_ref_grch38_20bp_XY.bed
    OUT=/storage/nmmsv/reference_TRF/hg38_ver${VER}_tmp.bed
fi
echo ">> Processing raw file: $TRF_IN"
TMP=/storage/nmmsv/reference_TRF/tmp/
DISEASE_BED=/storage/nmmsv/GangSTR/experimental/supplement_bed/${REF}_disease.txt
CODIS_BED=/storage/nmmsv/GangSTR/experimental/supplement_bed/${REF}_codis.bed
REMOVE_TMP=false

# Minimal trimming step
echo ">> Step1: Minimal trim. Results in: ${TMP}/${REF}_mintrim.bed"
python minimal_trim.py $TRF_IN ${TMP}/${REF}_mintrim.bed
echo ">> Step1: Done!"

# Removing hompolymer run motifs
echo ">> Step2: Remove homopolymer run motifs. Results in: ${TMP}/${REF}_mintrim_nohom.bed"
cat ${TMP}/${REF}_mintrim.bed | awk '$4 != 1' > ${TMP}/${REF}_mintrim_nohom.bed
echo ">> Step2: Done!"

# Removing bundles
echo ">> Step3: Remove region bundles. Results in: ${TMP}/${REF}_mintrim_nohom_nobundle.bed"
THRESH=50
python remove_bundles.py ${TMP}/${REF}_mintrim_nohom.bed ${TMP}/${REF}_mintrim_nohom_nobundle.bed $THRESH
echo ">> Step3: Done!"

# Removing messy regions
echo ">> Step4: Remove messy regions. Results in: ${TMP}/${REF}_mintrim_nohom_nobundle_clean.bed"
python remove_messy.py ${TMP}/${REF}_mintrim_nohom_nobundle.bed ${TMP}/${REF}_mintrim_nohom_nobundle_clean.bed
echo ">> Step4: Done!"

# Supplement with disease loci
echo ">> Step5: Supplement with disease loci. Results in: ${TMP}/${REF}_mintrim_nohom_nobundle_clean_supp.bed"
python supplement_ref.py ${TMP}/${REF}_mintrim_nohom_nobundle_clean.bed ${TMP}/${REF}_mintrim_nohom_nobundle_clean_supp.bed $DISEASE_BED,$CODIS_BED
echo ">> Step5: Done!"

# Removing extra columns
echo ">> Step6: Removing extra columns. Results in: $OUT"
cat ${TMP}/${REF}_mintrim_nohom_nobundle_clean_supp.bed | cut -f1-5 > $OUT
echo ">> Step6: Done!"

if [ "$REMOVE_TMP" = true ]; then
    echo ">> Removing temp files in $TMP"
    rm ${TMP}/${REF}_mintrim.bed
    rm ${TMP}/${REF}_mintrim_nohom.bed
    rm ${TMP}/${REF}_mintrim_nohom_nobundle.bed
    rm ${TMP}/${REF}_mintrim_nohom_nobundle_clean.bed
    rm ${TMP}/${REF}_mintrim_nohom_nobundle_clean_supp.bed
fi

echo
echo ">> All done!"



