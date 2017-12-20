#!/bin/bash

CHROMFA=$1
OUTFILE=$2
TRF=$3
NUMPROC=$4

BASE=$(basename -- "$0")
usage()
{
    BASE=$(basename -- "$0")
    echo "Generate a GangSTR reference file
Usage:
    $BASE ${CHROMFA} ${OUTFILE} ${TRF} ${NUMPROC}
    ${CHROMFA} is a directory of fasta files (chr1.fa, chr2.fa, ...chr22.fa)
    ${OUTFILE} path to output file
    ${TRF} points to the Tandem Repeat Finder binary
    ${NUMPROC} gives number of processors to use at once
"
    exit 1
}

die()
{
    echo "[$BASE] error: $@" >&2
    exit 1
}

echo "[$BASE] Running with options:"
echo "[$BASE]   CHROMFA=${CHROMFA}"
echo "[$BASE]   OUTFILE=${OUTFILE}"
echo "[$BASE]   TRF=${TRF}"
echo "[$BASE]   NUMPROC=${NUMPROC}"


# Check inputs
test -z ${CHROMFA} && usage
test -d ${CHROMFA} || die "${CHROMFA} is not a directory"
test -z ${OUTFILE} && usage
test -z ${TRF} && usage
if ! type ${TRF} > /dev/null 2>&1; then
    die "${TRF} not executable"
fi

matchscore=2
mismatchscore=10
indelscore=7
maxperiod=10 # Largest repeat unit
pm=80
pi=10
minscore=24 # Require at least 12 bp perfect matching
maxlen=1000

# Put temp outputs here
tmpdir=$(mktemp -d)
cd ${tmpdir}

echo "Running in ${tmpdir}"

for chrom in $(seq 1 22)
do
    echo "${TRF} \
	${CHROMFA}/chr${chrom}.fa \
	${matchscore} ${mismatchscore} ${indelscore} \
	${pm} ${pi} \
	${minscore} \
	${maxperiod} \
	-d -h -l 1"
done | xargs -P${NUMPROC} -I% -n 1 sh -c "%"

# TODO implement here any other filters we would like to have on the reference
for chrom in $(seq 1 22)
do
    cat ${tmpdir}/chr${chrom}.fa.${matchscore}.${mismatchscore}.${indelscore}.${pm}.${pi}.${minscore}.${maxperiod}.dat | \
	awk -F' ' '(NF==15)' | \
	awk -v "chrom=$chrom" -F' ' '{print "chr" chrom "\t" $1 "\t" $2 "\t" $3 "\t" $14 "\t" $15}' | \
	awk -v "maxlen=${maxlen}" '(($3-$2) <= maxlen)'
done > ${OUTFILE}
