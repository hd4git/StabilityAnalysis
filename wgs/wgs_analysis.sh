### Run command: ./wgs_br1d2.sh -sample <> -des <> -loc <> -res <> > >(tee logs/wgs_br1d2.err) 2> >(tee logs/wgs_br1d2.log >&2)
export JAVA_HOME=/usr/local/bioinf/java/jre1.8.0_91/
export PATH=$PATH:/usr/local/bioinf/java/jre1.8.0_91/bin

### Reference
ref="/data/borth/hdhiman/ref/picr_refseq/picr.fa"
gf="/data/borth/hdhiman/ref/picr_refseq/GCF_003668045.1_CriGri-PICR_genomic.gff"

### Directories
sample=$2
des=$4
results="/data/borth/hdhiman/janssen/results/wgs"
out=$results/$sample\_res
tmp="$out/tmp"
trimmed="/data/borth/hdhiman/samples/janssen/raw/genome/trimmed"

### Tools
fastqc="/usr/local/bioinf/fastqc/fastqc-0.11.5/fastqc"
trimmomatic="/usr/local/bioinf/trimmomatic/trimmomatic-0.36/trimmomatic-0.36.jar"
bwa=/usr/local/bioinf/bwa/bwa-0.7.15/bwa
java=/usr/local/bioinf/java/jre1.8.0_91/bin/java
picard=/usr/local/bioinf/picard/picard-tools-2.3.0
samtools=/usr/local/bioinf/samtools/samtools-1.9/bin/samtools
varscan=/usr/local/bioinf/varscan/varscan-2.4.2
gatk="/usr/local/bioinf/java/jre1.8.0_91/bin/java -jar /usr/local/bioinf/gatk/gatk-4.1.2.0/gatk-package-4.1.2.0-local.jar"
coovar="/data/borth/hdhiman/toolkit/CooVar-master/"

### Files
sample_trim=$trimmed/$sample"_R1_001_paired.fastq.gz "$trimmed/$sample"_R2_001_paired.fastq.gz"
sam_out=$out/wgs_$sample\.sam
bam_out=$out/wgs_$sample\.bam
sorted_bam_out=$out/wgs_$sample\_sorted.bam
rmdup_bam_out=$out/wgs_$sample\_sorted-dup.bam
final_rg_out=$out/wgs_$sample\_sorted-dup+rg.bam
final_rg_out_recal=$out/recal_data.table
gatk_out=$out/wgs_$sample\_var_GATK.txt
gatk_out_recal=$out/wgs_$sample\_var_GATK_recal.txt
doc=$out/wgs_$sample\-DOC
gatk_out=$out/wgs_$sample\_var_GATK.txt
realignment_targets_list=$out/realignment_targets_$sample\.list
realigned_reads_bam=$out/realigned_reads_$sample\.bam
raw_variants_vcf=$out/raw_variants_$sample\.vcf
raw_snps_vcf=$out/raw_snps_$sample\.vcf
raw_indels_vcf=$out/raw_indels_$sample\.vcf
filtered_snps_vcf=$out/filtered_snps_$sample\.vcf
filtered_indels_vcf=$out/filtered_indels_$sample\.vcf
recal_data_table=$out/recal_data_$sample\.table
post_recal_data_table=$out/post_recal_data_$sample\.table
recalibration_plots_pdf=$out/recalibration_plots_$sample\.pdf
recal_reads_bam=$out/recal_reads_$sample\.bam
raw_variants_recal_vcf=$out/raw_variants_recal_$sample\.vcf
raw_snps_recal_vcf=$out/raw_snps_recal_$sample\.vcf
raw_indels_recal_vcf=$out/raw_indels_recal_$sample\.vcf
filtered_indels_recal_vcf=$out/filtered_indels_recal_$sample\.vcf
filtered_snps_final_vcf=$out/filtered_snps_recal_$sample\.vcf
filtered_snv_final_vcf=$out/filtered_snv_recal_$sample\.vcf
filtered_snv_final_biallelic_vcf=$out/filtered_snv_recal_biallelic_$sample\.vcf
flagstat_ori=$out/wgs_$sample\_sorted_flagstat.txt
flagstat_final=$out/wgs_$sample\_sorted-dup+rg_flagstat.txt
filtered_var_final_vcf="$out/filtered_var_recal_$sample\.vcf.gz"

id=$sample
RGID=$sample
RGLB=$sample\_LB
RGPL="Illumina"
RGPU="HiSeq" 
RGSM=$sample
RGCN="Janssen" 
RGDS=$des

/usr/local/bioinf/bwa/bwa-0.7.15/bwa index -p picr -a bwtsw GCF_003668045.1_CriGri-PICR_genomic.fa
/usr/local/bioinf/samtools/samtools-1.9/bin/samtools faidx picr.fa

mkdir -p $tmp;
cd $out

#### Step 1. Alignment ####
$bwa mem -t 10 -M /data/borth/hdhiman/ref/picr_refseq/picr $sample_trim | $samtools view -bS - | $samtools sort - -o $sorted_bam_out;
    
#### Step 1''. Check alignment statistics for main alignment map ####
$samtools flagstat $sorted_bam_out > $flagstat_ori;
     
#### Step 2. Mark duplicates ####
$java -jar $picard/picard.jar  MarkDuplicates INPUT=$sorted_bam_out OUTPUT=$rmdup_bam_out METRICS_FILE=$out/metrics.txt REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE TMP_DIR=$tmp MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000;
      
#### Step 2'. Remove unprocessed alignment map ####
rm -rf $sorted_bam_out;
       
#### Step 3. Add read groups ####
$java -jar $picard/picard.jar AddOrReplaceReadGroups INPUT=$rmdup_bam_out OUTPUT=$final_rg_out RGID=$RGID RGLB=$RGLB RGPL=$RGPL RGPU=$RGPU RGSM=$RGSM RGCN=$RGCN RGDS=$RGDS CREATE_INDEX=true;
       
#### Step 3'. Remove alignment map without read groups ####
rm -rf $rmdup_bam_out;
        
#### Step 3". Check alignment statistics for alignment map after removing duplicates ####
$samtools flagstat $final_rg_out > $flagstat_final;
        
#### Step 4. Variant calling with GATK sec 10, scc 30 ####
#Call Variants
$gatk HaplotypeCaller -output-mode EMIT_ALL_CONFIDENT_SITES -stand-call-conf 30 -R $ref -I $final_rg_out -O $raw_variants_vcf
#Extract SNPs & Indels
$gatk SelectVariants -R $ref -V $raw_variants_vcf -select-type SNP -O $raw_snps_vcf
$gatk SelectVariants -R $ref -V $raw_variants_vcf -select-type INDEL -O $raw_indels_vcf
#Filter SNPs
$gatk VariantFiltration -R $ref -V $raw_snps_vcf --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0' --filter-name "basic_snp_filter" -O $filtered_snps_vcf
#Filter Indels
$gatk VariantFiltration -R $ref -V $raw_indels_vcf --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filter-name "basic_indel_filter" -O $filtered_indels_vcf
#Base Quality Score Recalibration (BQSR) 
$gatk BaseRecalibrator -R $ref -I $final_rg_out --known-sites $filtered_snps_vcf --known-sites $filtered_indels_vcf -O $recal_data_table
#Apply BQSR
$gatk ApplyBQSR -R $ref -I $final_rg_out --bqsr-recal-file $recal_data_table -O $recal_reads_bam
#Call Variants
$gatk HaplotypeCaller --output-mode EMIT_ALL_CONFIDENT_SITES -stand-call-conf 30 -R $ref -I $recal_reads_bam -O $raw_variants_recal_vcf
#Extract SNPs & Indels
$gatk SelectVariants -R $ref -V $raw_variants_recal_vcf -select-type SNP -O $raw_snps_recal_vcf
$gatk SelectVariants -R $ref -V $raw_variants_recal_vcf -select-type INDEL -O $raw_indels_recal_vcf
#Filter SNPs
$gatk VariantFiltration -R $ref -V $raw_snps_recal_vcf --filter-expression 'QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0' --filter-name "basic_snp_filter" -O $filtered_snps_final_vcf
# #Filter Indels
$gatk VariantFiltration -R $ref -V $raw_indels_recal_vcf --filter-expression 'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || SOR > 10.0' --filter-name "basic_indel_filter" -O $filtered_indels_recal_vcf

$java -jar /usr/local/bioinf/picard/picard-tools-2.3.0/picard.jar MergeVcfs I=$filtered_snps_final_vcf I=$filtered_indels_recal_vcf O=$filtered_var_final_vcf

#### Step 5. Estimate Depth of coverage using GATK ####
/usr/local/bioinf/qualimap/qualimap_v2.2.1/qualimap bamqc -bam $final_rg_out -outdir $out/cvg -outfile report_$sample\.pdf -outformat "PDF"


