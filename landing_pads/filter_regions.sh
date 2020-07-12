###################################################################################################
##############################	Post processing	###################################################
###################################################################################################

# check 20kb bins for distance from HighExpressionPeak, closest enhancer, promoter, ActiveTranscription mark, low variability region, repressed marks, quiescent regions
# check bins for presence in intron, exon or intergenic regions

while read f; 
do 
echo $f; 
perl scripts2/filter_regions.pl -in $f -win 20000 -version v1 &
done < src2/ref_split/list.txt

while read f; 
do 
echo $f; 
perl scripts2/filter_regions.pl -in $f -win 5000 -version v2 &
done < src2/ref_split/list.txt

# merge all result files
# add chr number
awk 'NR==FNR {a[$1]=$3; next} ($1 in a) {OFS="\t"; print $0,a[$1]}' \
	src2/picr.chromosome_assignment.txt \
	<(cat results2/filter_regions_v1_ref_split*) | \
	sort -k1,1 -k2,2n \
> results2/filter_regions_v1_all.txt

# replace "" with "_"
sed -i 's/\t\t/\t_\t/g' results2/filter_regions_v1_all.txt
sed -i 's/\t\t/\t_\t/g' results2/filter_regions_v1_all.txt
sed -i 's/\t\t/\t_\t/g' results2/filter_regions_v1_all.txt

# #filter bins for presence in lncRNA
# filter bins with no enhancer, promoter, active transcription mark or high expression peak in the scaffold
# filter bins inside exons, within high expression peak and with >250kb distance from high expression peak
# filter bins with <100kb distance from high variability region and <100kb distance from repressive marks 
# filter bins with lesser distance from repressor than all the active marks

###################################################################################################
###################################################################################################
###################################################################################################
# cat /var/www/html/jbrowse/data/sv_bed/df.*.sv.bed | cut -f1-3 | sort -k1,1 -k2,2n | uniq > /var/www/html/jbrowse/data/sv_bed/sv_all.bed
# awk '{if($2<0) {start=0} else {start=$2} print $1"\t"start"\t"$3}' /var/www/html/jbrowse/data/sv_bed/sv_all.bed > /var/www/html/jbrowse/data/sv_bed/sv_all_c.bed 
# ../toolkit/bedtools2-master/bin/bedtools merge -d 200 -i /var/www/html/jbrowse/data/sv_bed/sv_all_c.bed > /var/www/html/jbrowse/data/sv_bed/sv_all_merged.bed 
# Chr	Start	Stop	HighExpressionPeak	HighVariabilityPeak	Quiescent	Repressed	Enhancer	Promoter	ActiveTranscription	GeneRegion
intersectBed -v -a \
<(awk '{if($11 !~ /exon/ && $4>0 && $4<=250000 && $4!~ /_/ && \
$5>=3000 &&\
$7>$8 && $7>$9 && $7>$10 && \
$8!~ /_/ && $9!~ /_/ && $10!~ /_/) \
print $0}' \
results2/filter_regions_v1_all.txt) \
-b src2/sv_all_merged.bed \
> results2/filter_regions_v1_all_filtered.txt


# merge closeby bins 
bedtools merge -d 2001 -i results2/filter_regions_v1_all_filtered.txt |\
sort -k1,1 -k2,2n | uniq \
> results2/filter_regions_v1_all_parsed.txt

# get information for merged bins
perl scripts2/filter_regions2.pl -in filter_regions_v1_all_parsed

# add chromosome number 
# sort 

awk 'NR==FNR {a[$1]=$3; next} ($1 in a) \
	{OFS="\t"; print $0"\t"a[$1]"\t"$3-$2}' \
	src2/picr.chromosome_assignment.txt \
	results2/filter_regions_v1_all_parsed_filtered.txt |
	sort -k5,5hr -k7,7hr -k4,4h -k10,10h -k13,13h -k12,12h |  awk '{print $0"\t"FNR}' > results2/filter_regions_v1_all_parsed_filtered_sorted.txt


awk 'NR==FNR {a[$1]=$3; next} ($1 in a) \
	{OFS="\t"; print $0"\t"a[$1]"\t"$3-$2}' \
	src2/picr.chromosome_assignment.txt \
	results2/filter_regions_v1_all_parsed_filtered.txt \
	> results2/filter_regions_v1_all_parsed_filtered2.txt

intersectBed -a results2/filter_regions_v1_all_parsed_filtered_sorted.txt -b src2/unfavorable_ncbi.txt -wa -wb > results2/intersection.unfav.txt
intersectBed -a results2/filter_regions_v1_all_parsed_filtered_sorted.txt -b src2/favorable_ncbi.txt -wa -wb > results2/intersection.fav.txt

#################################################################################################################################
# filter bins with no active/inactive states in the scaffold
# filter bins for presence either inside exons or with high expression peak >250kb distance 
# filter bins for presence of high variability region <50kb distance and <50kb distance from repressive marks 
# filter bins with lesser distance from repressor than all the active marks

# merge all result files
# add chr number
awk 'NR==FNR {a[$1]=$2; next} ($1 in a) {OFS="\t"; print $0,a[$1]}' \
src2/picr.chromosome_assignment.txt \
<(cat results2/filter_regions_v2_ref_split*) | \
sort -k1,1 -k2,2n \
> results2/filter_regions_v2_all.txt

# replace "" with "_"
sed -i 's/\t\t/\t_\t/g' results2/filter_regions_v2_all.txt
sed -i 's/\t\t/\t_\t/g' results2/filter_regions_v2_all.txt
sed -i 's/\t\t/\t_\t/g' results2/filter_regions_v2_all.txt

intersectBed -v -a <(awk '{if($5<=500 || ($7<=$8 && $7<=$9 && $7<=$10 && $4>=700000 && $10>500000)) print $0}' \
<(grep -v -P "\t_\t" results2/filter_regions_v2_all.txt) ) -b src2/sv_all_merged.bed \
> results2/filter_regions_v2_all_filtered.txt

# merge closeby bins 
../../../toolkit/bedtools2-master/bin/bedtools merge -d 10001 -i results2/filter_regions_v2_all_filtered.txt |\
sort -k1,1 -k2,2n | uniq \
> results2/filter_regions_v2_all_parsed.txt

# get information for merged bins
perl scripts2/filter_regions2.pl -in filter_regions_v2_all_parsed

awk 'NR==FNR {a[$1]=$3; next} ($1 in a) \
	{OFS="\t"; print $0"\t"a[$1]"\t"$3-$2}' \
	src2/picr.chromosome_assignment.txt \
	<(awk '{if($5<=500 || ($7<=$8 && $7<=$9 && $7<=$10 && $4>=700000 && $10>500000)) print $0}' results2/filter_regions_v2_all_parsed_filtered.txt) \
> results2/filter_regions_v2_all_parsed_filtered2.txt

### further analysis with filter_regions.R
