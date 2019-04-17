#!/bin/bash

REF=$1
TRF_IN=$2
VER=$3
OUT_DIR=$4

BASE=$(basename -- "$0")
usage()
{
    BASE=$(basename -- "$0")
    echo "Trim a GangSTR reference file
Usage:
    $BASE {REF} {TRF_IN} {TRF} {NUMPROC}
    {REF} is a for the reference (hg38, hg18, mm10, etc.)
    {TRF_IN} path to results from TRF
    {VER} version number
    {OUT_DIR} path to output directory.
"
    exit 1
}

die()
{
    echo "[$BASE] error: $@" >&2
    exit 1
}

echo "[$BASE] Running with options:"
echo "[$BASE]   REF=${REF}"
echo "[$BASE]   TRF_IN=${TRF_IN}"
echo "[$BASE]   VER=${VER}"
echo "[$BASE]   OUT_DIR=${OUT_DIR}"



# Check inputs
test -e ${TRF_IN} || die "cannot access ${TRF_IN}: No such file or directory"
test -d ${OUT_DIR} || die "cannot access ${OUT_DIR}: No such file or directory"


REMOVE_TMP=false
SUPP_FLAG=false
DISEASE_BED=/storage/nmmsv/GangSTR/experimental/supplement_bed/${REF}_disease.txt
CODIS_BED=/storage/nmmsv/GangSTR/experimental/supplement_bed/${REF}_codis.bed

#if [ "$REF" = "hg19" ]; then
#    TRF_IN=/storage/nmmsv/reference_TRF/raw/raw_ref_hg19_20bp_XY.bed
#    OUT=/storage/nmmsv/reference_TRF/hg19_ver${VER}_tmp.bed
#elif [ "$REF" = "hg38" ]; then
#    TRF_IN=/storage/nmmsv/reference_TRF/raw/raw_ref_grch38_20bp_XY.bed
#    OUT=/storage/nmmsv/reference_TRF/hg38_ver${VER}_tmp.bed
#fi

echo ">> Processing raw file: $TRF_IN"
TMP=$OUT_DIR/tmp/
mkdir -p $TMP
OUT=${OUT_DIR}/${REF}_ver${VER}.bed

# Minimal trimming step
echo ">> Step1: Minimal trim. Results in: ${TMP}/${REF}_mintrim.bed"
python minimal_trim.py $TRF_IN ${TMP}/${REF}_mintrim.bed
echo ">> Step1: Done!"

# Removing hompolymer run motifs

echo ">> Step2: Remove homopolymer run motifs. Results in: ${TMP}/${REF}_mintrim_nohom.bed"
cat ${TMP}/${REF}_mintrim.bed | awk '$4 != 1' > ${TMP}/${REF}_mintrim_nohom.bed
#cat ${TMP}/${REF}_mintrim.bed > ${TMP}/${REF}_mintrim_nohom.bed
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
if [ "$SUPP_FLAG" = true ]; then
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
	rm -r ${TMP}
    fi
else
    echo ">> Step5: SKIPPED!"
    # Removing extra columns
    echo ">> Step6: Removing extra columns. Results in: $OUT"
    cat ${TMP}/${REF}_mintrim_nohom_nobundle_clean.bed | cut -f1-5 > $OUT
    echo ">> Step6: Done!"

    if [ "$REMOVE_TMP" = true ]; then
	echo ">> Removing temp files in $TMP"
	rm ${TMP}/${REF}_mintrim.bed
	rm ${TMP}/${REF}_mintrim_nohom.bed
	rm ${TMP}/${REF}_mintrim_nohom_nobundle.bed
	rm ${TMP}/${REF}_mintrim_nohom_nobundle_clean.bed
	rm -r ${TMP}
    fi
fi



echo
echo ">> All done!"



