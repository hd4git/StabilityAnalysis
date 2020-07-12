### Run command: /data/borth/hdhiman/scripts/janssen/wgs/./trim.sh  <sample ID> > >(tee /data/borth/hdhiman/scripts/janssen/wgs/trim_<ID>.err) 2> >(tee /data/borth/hdhiman/scripts/janssen/wgs/trim_<ID>.log >&2)

bowtie2="/usr/local/bioinf/bowtie2/bowtie2-2.2.9/bin/bowtie2"
f=$1
echo $f
mkdir filtered/$f

### Map to confidential vector sequence and use unmapped reads for further analysis
$bowtie2 -p 10 \
		--no-unal \
		--un-conc-gz filtered/$f/ \
		--al-conc-gz filtered/$f/ \
		-x ../vector_sequence/VectorA \
		-1 $f\_R1_001.fastq.gz \
		-2 $f\_R2_001.fastq.gz | \
    	samtools sort - | \
    	samtools view -bS > vector_bam/$f\_R1_mapped+unmapped.bam

### Rename files
mv filtered/$f//un-conc-mate.1 filtered/$f/$f\_R1_001_filtered.fastq.gz
mv filtered/$f/un-conc-mate.2 filtered/$f/$f\_R2_001_filtered.fastq.gz
mv filtered/$f/al-conc-mate.1 filtered/$f/$f\_R1_001_aligned.fastq.gz
mv filtered/$f/al-conc-mate.2 filtered/$f/$f\_R2_001_aligned.fastq.gz

### QC raw data
/usr/local/bioinf/fastqc/fastqc-0.11.8/./fastqc -o fastqc_filtered filtered/$f/$f\_R1_001_filtered.fastq.gz filtered/$f/$f\_R2_001_filtered.fastq.gz

### Trim filtered files
/usr/local/bioinf/java/jre1.8.0_91/bin/java -jar /usr/local/bioinf/trimmomatic/trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 \
    filtered/$f/$f\_R1_001_filtered.fastq.gz \
    filtered/$f/$f\_R2_001_filtered.fastq.gz \
    trimmed/$f\_R1_001_paired.fastq.gz \
    trimmed/$f\_R1_001_unpaired.fastq.gz \
    trimmed/$f\_R2_001_paired.fastq.gz \
    trimmed/$f\_R2_001_unpaired.fastq.gz \
    ILLUMINACLIP:/usr/local/bioinf/trimmomatic/trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:25

### QC filtered data
/usr/local/bioinf/fastqc/fastqc-0.11.8/./fastqc -o fastqc_trimmed trimmed/$f\_R1_001_paired.fastq.gz trimmed/$f\_R2_001_paired.fastq.gz


