source /data/borth/hdhiman/scripts/janssen/var_list-rnaseq.conf
out="$results/G9_res"
tmp="$out/tmp"
name="G9"

mkdir -p $tmp;

## Replicate 1
$hisat2 -p 10 -x $index --rg-id $name --rg SM:Replicate1 $sample_G9_fq1 | $samtools view -bS - | $samtools sort - -o $out/$name\_rep1.bam 
/usr/local/bioinf/python/python-2.7.14/bin/htseq-count -f bam -i gene -s no $out/$name\_rep1.bam $gff > $out/$name\_rep1.counts

## Replicate 2
$hisat2 -p 10 -x $index $sample_G9_fq2 --rg-id $name --rg SM:Replicate2  | $samtools view -bS - | $samtools sort - -o $out/$name\_rep2.bam
/usr/local/bioinf/python/python-2.7.14/bin/htseq-count -f bam -i gene -s no $out/$name\_rep2.bam $gff > $out/$name\_rep2.counts

## Replicate 3
$hisat2 -p 10 -x $index --rg-id $name --rg SM:Replicate3 $sample_G9_fq3 | $samtools view -bS - | $samtools sort - -o $out/$name\_rep3.bam
/usr/local/bioinf/python/python-2.7.14/bin/htseq-count -f bam -i gene -s no $out/$name\_rep3.bam $gff > $out/$name\_rep3.counts

### Run command./rnaseq_G9.sh > >(tee logs/rnaseq_G9.err) 2> >(tee logs/rnaseq_G9.log >&2)