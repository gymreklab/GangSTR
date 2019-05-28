sam=$1
out=$2
#chrom=$3
#start=$4
#end=$5


vicin=100000

samtools view $sam | cut -f3,4 | uniq | awk '{print $0"\t"$2}' > $out

# Filtering regions in vicinity:
#samtools view $sam | awk -v awkchrom="$chrom" -v awkstart=$start -v awkend=$end -v awkvicin=$vicin\
# '($3 != awkchrom) || (($3 == awkchrom) && (($4 < awkstart - awkvicin) || ($4 > awkend + awkvicin))) {print}' | \
#cut -f3,4 | uniq | awk '{print $0"\t"$2}' > $out


