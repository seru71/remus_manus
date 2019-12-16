# link in REMUS/data/REG for REG in {promoter, enhancers, chromatin}:
# ln -s REMUS/data/promoters .

# generate all tissue BEDs
for reg in enhancers promoters chromatin; do
	for dir in $reg/*/[hG]*; do
		echo Processing ${dir}
		mkdir -p alltissues/${dir}
		zcat $dir/*.gz | sort -k 1V -k 2n | bedtools merge -i - > alltissues/${dir}/all.bed
	done
done


for f in alltissues/*/*/*/all.bed; do
	echo -n $f | sed 's/\//\t/g'
        awk '{sum+=$3-$2}END{print "\t"sum"\t"NR}' $f
done | cut -f2- > alltissues.tsv

function get_stats {
	echo -n $1 | sed 's/\//\t/g' | sed 's/\(_[0-9]*\)_/\1\t/' | sed 's/.bed.gz//'
        zcat $1 | awk '{sum+=$3-$2}END{print "\t"sum"\t"NR}'
}

for f in chromatin/*/[Gh]*/*.gz enhancers/*/[Gh]*/*.gz promoters/*/[Gh]*/*.gz; do 
	get_stats $f
done > reg_region_stats.tsv


