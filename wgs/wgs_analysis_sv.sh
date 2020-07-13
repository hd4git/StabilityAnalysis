### Run command: ./wgs_analysis_sv.sh -sample <sample> > >(tee logs/wgs_analysis_sv_<sample>.err) 2> >(tee logs/wgs_analysis_sv_<sample>.log >&2)

export JAVA_HOME=/usr/local/bioinf/java/latest/
export PATH=$PATH:/usr/local/bioinf/java/latest/bin
java=/usr/local/bioinf/java/latest/bin/java
delly="/usr/local/bioinf/delly/delly-0.7.9"
kentUtils="/usr/local/bioinf/kentUtils/src/kentUtils-302.1.0/samtabix"
samtools="/usr/local/bioinf/samtools/samtools-1.9/bin/samtools"
bedtools="/usr/local/bioinf/bedtools/bedtools-2.26.0/bin/bedtools"
bwa=/usr/local/bioinf/bwa/bwa-0.7.15/bwa
picard=/usr/local/bioinf/picard/picard-tools-2.3.0
samtools=/usr/local/bioinf/samtools/samtools-1.9/bin/samtools
varscan=/usr/local/bioinf/varscan/varscan-2.4.2
gatk="/data/borth/hdhiman/toolkit/gatk-4.1.1.0/gatk"
lumpy="/usr/local/bioinf/lumpy"

### Directories
sample=$2
results="/data/borth/hdhiman/janssen/results/wgs"
out=$results/$sample\_res
tmp="$out/tmp"

ref=/data/borth/hdhiman/ref/picr_refseq/picr.fa
mito_shifted="/data/borth/hdhiman/ref/c_griseus-mito_shifted.fasta"

final_rg_out=$out/wgs_$sample\_sorted-dup+rg.bam
recal_reads_bam=$out/recal_reads_$sample\.bam
raw_variants_recal_vcf=$out/raw_variants_recal_$sample\.vcf
raw_snps_recal_vcf=$out/raw_snps_recal_$sample\.vcf
raw_indels_recal_vcf=$out/raw_indels_recal_$sample\.vcf
filtered_indels_recal_vcf=$out/filtered_indels_recal_$sample\.vcf
filtered_snps_final_vcf=$out/filtered_snps_recal_$sample\.vcf


###########################################################################
#####################	Delly 	###########################################
###########################################################################
$delly/delly-0.7.9 call -r 20 -g $ref -o $out/wgs_$sample\_sv-n.bcf $recal_reads_bam
$kentUtils/bcftools/bcftools view -i 'MIN(INFO/PE>4) & MIN(INFO/MAPQ>19)' $out/wgs_$sample\_sv-n.bcf > $out/wgs_$sample\_sv-n.vcf
$kentUtils/bgzip $out/wgs_$sample\_sv-n.vcf
$kentUtils/tabix -p vcf $out/wgs_$sample\_sv-n.vcf.gz

###########################################################################
#####################	Lumpy 	###########################################
###########################################################################
# Extract the discordant paired-end alignments.
$samtools view -b -F 1294 $recal_reads_bam > $out/recal_reads_$sample\.discordants.unsorted.bam

# Extract the split-read alignments
$samtools view -h $recal_reads_bam \
   | $lumpy/lumpy_src/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
   | $samtools view -Sb - \
   > $out/recal_reads_$sample\.splitters.unsorted.bam

# Sort both alignments
$samtools sort  -o $out/recal_reads_$sample\.discordants.sorted.bam $out/recal_reads_$sample\.discordants.unsorted.bam
$samtools sort -o  $out/recal_reads_$sample\.splitters.sorted.bam $out/recal_reads_$sample\.splitters.unsorted.bam

#Run LUMPY Express on a single sample with pre-extracted splitters and discordants
$lumpy/lumpyexpress \
    -B $recal_reads_bam \
    -S $out/recal_reads_$sample\.splitters.sorted.bam \
    -D $out/recal_reads_$sample\.discordants.sorted.bam \
    -o $out/wgs_$sample\_sv-lumpy.vcf

###########################################################################
#####################	Manta 	###########################################
###########################################################################

python /usr/local/bioinf/manta/manta-1.6.0_prog/bin/configManta.py --bam recal_reads_$sample\.bam --referenceFasta /data/borth/hdhiman/ref/picr_refseq/picr.fa  --runDir $out/manta_analysis
$out/manta_analysis/runWorkflow.py -m local -j 8
