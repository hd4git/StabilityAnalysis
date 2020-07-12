while read f;
do 
echo $f
#/data/borth/hdhiman/toolkit/bbmap/./bbsplit.sh  in1=unfiltered/$f\-3_R1_001.fastq.gz in2=unfiltered/$f\-3_R2_001.fastq.gz ref=../vector_sequence/VectorA.fa basename=$f\-3_out_%.fq outu1=filtered/$f\-3_R1_001_filtered.fastq outu2=filtered/$f\-3_R2_001_filtered.fastq
/usr/local/bioinf/bowtie2/bowtie2-2.2.9/bin/bowtie2 -p 10 --no-unal --un-conc-gz filtered3/ --al-conc-gz filtered3/ -x ../vector_sequence/VectorA -1 unfiltered/$f\-3_R1_001.fastq.gz -2 unfiltered/$f\-3_R2_001.fastq.gz | samtools sort - | samtools view -bS > $f\-3_R1_mapped+unmapped.bam
mv filtered3/un-conc-mate.1 filtered3/$f\-3_R1_001_filtered.fastq.gz
mv filtered3/un-conc-mate.2 filtered3/$f\-3_R2_001_filtered.fastq.gz
mv filtered3/al-conc-mate.1 filtered3/$f\-3_R1_001_aligned.fastq.gz
mv filtered3/al-conc-mate.2 filtered3/$f\-3_R2_001_aligned.fastq.gz
/usr/local/bioinf/fastqc/fastqc-0.11.8/./fastqc -o fastqc_filtered filtered3/$f\-3_R1_001_filtered.fastq.gz filtered3/$f\-3_R2_001_filtered.fastq.gz
/usr/local/bioinf/java/jre1.8.0_91/bin/java -jar /usr/local/bioinf/trimmomatic/trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 filtered3/$f\-3_R1_001_filtered.fastq.gz filtered3/$f\-3_R2_001_filtered.fastq.gz trimmed/$f\-3_R1_001_paired.fastq.gz trimmed/$f\-3_R1_001_unpaired.fastq.gz trimmed/$f\-3_R2_001_paired.fastq.gz trimmed/$f\-3_R2_001_unpaired.fastq.gz  ILLUMINACLIP:/usr/local/bioinf/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25
/usr/local/bioinf/fastqc/fastqc-0.11.8/./fastqc -o fastqc_trimmed trimmed/$f\-3_R1_001_paired.fastq.gz trimmed/$f\-3_R2_001_paired.fastq.gz
done < /data/borth/hdhiman/scripts/janssen/samples.txt

#/data/borth/hdhiman/scripts/janssen/./trim3.sh  > >(tee /data/borth/hdhiman/scripts/janssen/trim3.err) 2> >(tee /data/borth/hdhiman/scripts/janssen/trim3.log >&2)