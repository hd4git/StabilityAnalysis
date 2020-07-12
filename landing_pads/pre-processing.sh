## Bin the genome
## cd /data/borth/hdhiman/janssen/integration_sites/filter_regions/src2
perl ../scripts/genome_bin.pl -bin 1000 -index /data/borth/hdhiman/ref/picr_refseq/picr.fa.fai
perl ../scripts/genome_bin.pl -bin 2000 -index /data/borth/hdhiman/ref/picr_refseq/picr.fa.fai
perl ../scripts/genome_bin.pl -bin 3000 -index /data/borth/hdhiman/ref/picr_refseq/picr.fa.fai
perl ../scripts/genome_bin.pl -bin 5000 -index /data/borth/hdhiman/ref/picr_refseq/picr.fa.fai
perl ../scripts/genome_bin.pl -bin 8000 -index /data/borth/hdhiman/ref/picr_refseq/picr.fa.fai
perl ../scripts/genome_bin.pl -bin 10000 -index /data/borth/hdhiman/ref/picr_refseq/picr.fa.fai

################################################################################################
## Variant calls to bed
while read f; 
do 
	zgrep -v "#" /data/borth/hdhiman/janssen/results/wgs/$f\_res/filtered_var_recal_$f\.vcf.gz |\
	awk '{len=length($4); print $1"\t"$2"\t"$2+len}' \
	> variant_calls/filtered_var_recal_$f\.bed & \
done < list.txt

## Merge all
cat variant_calls/*bed | sort -k1,1 -k2,2n | uniq > variant_calls/all_small_variations.bed

## Get binned counts
while read f; 
do 
/usr/local/bioinf/bedtools/latest/intersectBed -a variant_calls/all_small_variations.bed -b genome_binned-$f\000.bed  -wa -wb  | cut -f4-6 | sort | uniq -c | awk '{OFS="\t"; print $2,$3,$4,$1}' | sort -k1,1 -k2,2n > variant_calls/bg/all_snp-$f\kb.bedGraph &
done < list_bin_size.txt

while read f; 
do 
/usr/local/bioinf/bedtools/latest/intersectBed -a variant_calls/filtered_var_recal_$f\.bed -b genome_binned-2000.bed  -wa -wb  | cut -f4-6 | sort | uniq -c | awk '{OFS="\t"; print $2,$3,$4,$1}' | sort -k1,1 -k2,2n > variant_calls/bg/small_variations_$f\-2kb.bedGraph &
done < list.txt


## Get M-Values for variant counts
#./Rscript mval_var.R

## Conver bedGraphToBigWig for all variant call for binned genome (M-values)
while read f; 
do 
/usr/local/bioinf/kentUtils/latest/bedGraphToBigWig variant_calls/bg/$f\_mval.bedGraph /data/borth/hdhiman/ref/picr_refseq/chrom.sizes variant_calls/bw/$f\_mval.bw &
done < /data/borth/hdhiman/janssen/integration_sites/filter_regions/src2/list_var_bg.txt


################################################################################################
## Index rna-seq bam files
while read f; do echo $f; samtools index $f\.bam & done < bam_list1.txt 
while read f; do echo $f; samtools index $f\.bam & done < bam_list2.txt 
while read f; do echo $f; samtools index $f\.bam & done < bam_list3.txt 
	
## bedGraph of expression profiles
while read f; do echo $f; bamCoverage -b $f\.bam  -o $f\.bedGraph -of bedgraph -p 10; done < bam_list1.txt &
while read f; do echo $f; bamCoverage -b $f\.bam  -o $f\.bedGraph -of bedgraph -p 10; done < bam_list2.txt &
while read f; do echo $f; bamCoverage -b $f\.bam  -o $f\.bedGraph -of bedgraph -p 10; done < bam_list3.txt &

	
## Get M-Values for expression levels
#./Rscript mval_exp.R

## Conver bedGraphToBigWig for expression levels for binned genome (M-values)
while read f; 
do 
/usr/local/bioinf/kentUtils/latest/bedGraphToBigWig exp/bg/$f\_mval.bedGraph /data/borth/hdhiman/ref/picr_refseq/chrom.sizes exp/bw/$f\_mval.bw &
done < /data/borth/hdhiman/janssen/integration_sites/filter_regions/src2/list_exp_bg.txt

## Expression peaks for landing pads 
while read f; 
do 
echo $f; 
/usr/local/bioinf/bedtools/latest/bedtools merge -d 1 -i \
	<(awk '{if($4>=20) print $0}' exp/bg/$f\_mval.bedGraph) > \
	exp/bg_gt20/$f\_mval_gt20.bed ; 
done < list_exp_bg.txt
/usr/local/bioinf/bedtools/latest/bedtools merge -d 1 -i <(cat exp/bg_gt20/*.bed | sort -k1,1 -k2,2n) > exp/bg_gt20/exp_gt20.bed

