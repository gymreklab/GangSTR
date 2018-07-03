#!/bin/bash
echo "Running GangSTR test:"
ref=path/to/hg38/refgenome.fa

for bam in alignment/*.bam
do
    name=$(basename $bam .sorted.bam)
    echo Sample: $name
    ../src/GangSTR \
	--bam $bam \
	--ref $ref \
	--regions HTT.bed \
	--out results/$name \
	--coverage 50 \
	--useofftarget \
	--output-readinfo \
	--verbose
done
