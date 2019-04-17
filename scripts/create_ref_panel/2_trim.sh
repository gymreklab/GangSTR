#!/bin/bash

usage()
{
    BASE=$(basename -- "$0")
    echo "Trim a GangSTR reference file
Usage:
    $BASE {TRF_IN} {OUT_PATH}
    {TRF_IN} path to results from TRF
    {OUT_PATH} path to output file.
"
    exit 1
}

die()
{
    echo "[$BASE] error: $@" >&2
    exit 1
}

if [ $# -eq 0 ]
  then
    usage
fi
if [ -z "$1" ]
  then
    usage
fi
if [ -z "$2" ]
  then
    usage
fi

TRF_IN=$1
OUT_PATH=$2
OUT_DIR=$(dirname "${OUT_PATH}")
OUT_FILE=$(basename -- "$OUT_PATH")
OUT_PREF="${OUT_FILE%.*}"






BASE=$(basename -- "$0")


# Check inputs
test -e ${TRF_IN} || die "cannot access ${TRF_IN}: No such file or directory"
test -d ${OUT_DIR} || die "cannot access ${OUT_DIR}: No such file or directory"


REMOVE_TMP=false
TMP=$OUT_DIR/tmp/
mkdir -p $TMP


echo "[$BASE] Running with options:"
echo "[$BASE]   TRF_IN=${TRF_IN}"
echo "[$BASE]   OUT_PATH=${OUT_PATH}"
echo "[$BASE]   TMP=${TMP}"





echo ">> Processing raw file: $TRF_IN"


# Minimal trimming step
echo ">> Step1: Minimal trim. Results in: ${TMP}/${OUT_PREF}_mintrim.bed"
python scripts/minimal_trim.py $TRF_IN ${TMP}/${OUT_PREF}_mintrim.bed
echo ">> Step1: Done!"

# Removing hompolymer run motifs

echo ">> Step2: Remove homopolymer run motifs. Results in: ${TMP}/${OUT_PREF}_mintrim_nohom.bed"
cat ${TMP}/${OUT_PREF}_mintrim.bed | awk '$4 != 1' > ${TMP}/${OUT_PREF}_mintrim_nohom.bed
#cat ${TMP}/${REF}_mintrim.bed > ${TMP}/${REF}_mintrim_nohom.bed
echo ">> Step2: Done!"


# Removing bundles
echo ">> Step3: Remove region bundles. Results in: ${TMP}/${OUT_PREF}_mintrim_nohom_nobundle.bed"
THRESH=50
python scripts/remove_bundles.py ${TMP}/${OUT_PREF}_mintrim_nohom.bed ${TMP}/${OUT_PREF}_mintrim_nohom_nobundle.bed $THRESH
echo ">> Step3: Done!"

# Removing messy regions
echo ">> Step4: Remove messy regions. Results in: ${TMP}/${OUT_PREF}_mintrim_nohom_nobundle_clean.bed"
python scripts/remove_messy.py ${TMP}/${OUT_PREF}_mintrim_nohom_nobundle.bed ${TMP}/${OUT_PREF}_mintrim_nohom_nobundle_clean.bed
echo ">> Step4: Done!"

# Removing extra columns
echo ">> Step5: Removing extra columns. Results in: $OUT_PATH"
cat ${TMP}/${OUT_PREF}_mintrim_nohom_nobundle_clean.bed | cut -f1-5 > $OUT_PATH
echo ">> Step5: Done!"


if [ "$REMOVE_TMP" = true ]; then
    echo ">> Removing temp files in $TMP"
    rm ${TMP}/${OUT_PREF}_mintrim.bed
    rm ${TMP}/${OUT_PREF}_mintrim_nohom.bed
    rm ${TMP}/${OUT_PREF}_mintrim_nohom_nobundle.bed
    rm ${TMP}/${OUT_PREF}_mintrim_nohom_nobundle_clean.bed
    rm -r ${TMP}
fi


echo
echo ">> All done!"

exit 1

# Supplement with disease loci
if [ "$SUPP_FLAG" = true ]; then
    echo ">> Step5: Supplement with disease loci. Results in: ${TMP}/${REF}_mintrim_nohom_nobundle_clean_supp.bed"
    python supplement_ref.py ${TMP}/${REF}_mintrim_nohom_nobundle_clean.bed ${TMP}/${REF}_mintrim_nohom_nobundle_clean_supp.bed $DISEASE_BED,$CODIS_BED
    echo ">> Step5: Done!"
    
    
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






