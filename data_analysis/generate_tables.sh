# in REMUS/ run

# generate all tissue BEDs
#for reg in enhancers promoters chromatin; do
#	for dir in data/$reg/*/[hG]*; do
#		echo Processing ${dir}
#		mkdir -p data_analysis/${dir}
#		zcat $dir/*.gz | sort -k 1V -k 2n | bedtools merge -i - > data_analysis/${dir}/all.bed
#	done
#done


cd data_analysis
mv data alltissues
for f in alltissues/*/*/*/all.bed; do
	echo -n $f | sed 's/\//\t/g'
        awk '{sum+=$3-$2}END{print "\t"sum"\t"NR}' $f
done | cut -f2- > alltissues.tsv

function get_stats {
	echo -n $1 | sed 's/\//\t/g' | sed 's/\(_[0-9]*\)_/\1\t/' | sed 's/.bed.gz//'
        zcat $1 | awk '{sum+=$3-$2}END{print "\t"sum"\t"NR}'
}

cd ..
for f in data/chromatin/*/[Gh]*/*.gz data/enhancers/*/[Gh]*/*.gz data/promoters/*/[Gh]*/*.gz; do 
	get_stats $f
done > data_analysis/reg_region_stats.tsv


