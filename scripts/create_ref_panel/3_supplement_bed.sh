
IN_BED=/PATH/TO/TRIMMED.bed
OUT_BED=/PATH/TO/TRIMMED_supp.bed
SUPP_BED1=/PATH/TO/DISEASE.bed
SUPP_BED2=/PATH/TO/CODIS.bed
SUPP_BED3=/PATH/TO/ADDITIONAL.bed

echo ">> Supplement with extra loci."
python scripts/supplement_ref.py $IN_BED $OUT_BED $SUPP_BED1,$SUPP_BED2,$SUPP_BED3
echo ">> Step"
