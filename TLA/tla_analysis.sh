while read f;
do
/usr/local/bioinf/qualimap/qualimap_v2.2.1/qualimap --java-mem-size=30G bamqc -bam $f\.bam -outdir $f\_cvg -outfile $f\_report.pdf -outformat "PDF"  
done < bam_files.txt



bam2fq="/usr/local/bioinf/bedtools/bedtools-2.28.0/bin/bamToFastq"
ref="/data/borth/hdhiman/ref/picr_refseq/picr.fa"
samtools="/usr/local/bioinf/samtools/samtools-1.9/bin/samtools"
bwa="/usr/local/bioinf/bwa/bwa-0.7.15/bwa"
while read i;
do
f=(${i//_L/ });
echo $f;
# $bam2fq -i bam/$i -fq fastq/$f\.fastq 
# gzip fastq/$f\.fastq 
# $bwa bwasw -t20 -b7 /data/borth/hdhiman/ref/picr_refseq/picr fastq/$f\.fastq.gz | \
# 	$samtools view -bS - | \
# 	$samtools sort - -o out/$f\.bam
# $samtools index out/$f\.bam
/usr/local/bioinf/bedtools/bedtools-2.28.0/bin/bedtools genomecov -bg -ibam out/$f\.bam > out/$f\.bedGraph
done < list_bam.txt


bam2fq="/usr/local/bioinf/bedtools/bedtools-2.28.0/bin/bamToFastq"
ref="/data/borth/hdhiman/ref/picr_refseq/picr.fa"
samtools="/usr/local/bioinf/samtools/samtools-1.9/bin/samtools"
bwa="/usr/local/bioinf/bwa/bwa-0.7.15/bwa"
while read i;
do
f=(${i//_L/ });
echo $f;
$bam2fq -i bam/$i -fq fastq/$f\.fastq 
gzip fastq/$f\.fastq 
$bwa bwasw -t20 -b7 /data/borth/hdhiman/ref/picr_refseq/picr fastq/$f\.fastq.gz | \
	$samtools view -bS - | \
	$samtools sort - -o out/$f\.bam
done < bam_list2.txt



samtools merge -f <sample1.bam> <sample1_R1.bam> <sample1_R2.bam>
igvtools count -w 25 -f mean,median <sample1.bam> <sample1.bam.tdf> <ref_genome.fa>



##########################################
require("stringr")
require("data.table")
require("GenomicRanges")
require("Gviz")
require("genomation")
require("GenomicFeatures")

files<-list.files(pattern="*bedGraph$")
names<-as.data.frame(str_split_fixed(files, "_S", 2))
info<-lapply(files, function(x) {fread(x, header=F, col.names = c('chromosome', 'start', 'end', 'value')) } )
names(info)<-names$V1
setwd("../TLA")
options(ucscChromosomeNames=FALSE) 


# gff<-gffToGRanges("/home/hdhiman/Data/work/hd/ref/picr/ncbi/GCF_003668045.1_CriGri-PICR_genomic.gff",
# 					filter = NULL, zero.based = FALSE, ensembl = FALSE)
# gff<-gff[which(gff$source != "RefSeq"), c(1,2,21)]
# gff.df<-as.data.frame(gff)

# gtrack <- GenomeAxisTrack()
# atrack <- GeneRegionTrack(gff, chromosome="NW_020822638.1", start=1245001, end=2426000,
# 	cex.group=1,
# 	name = "Genes", 
# 	showID = TRUE, 
# 	geneSymbol = TRUE)
# plotTracks(c(gtrack,atrack),  transcriptAnnotation = "symbol", stacking = "squish")
# dev.off()

sample.G9 <- info[30:33]
sample.5G10 <- info[13:16]
sample.1E3 <- info[5:8]
sample.1C11 <- info[1:4]
sample.6A6 <- info[17:20]
sample.6H1 <- info[21:24]


# makeGRangesObject <- function(Chrom = NULL, Start = NULL, End = NULL, 
#     Strand = NULL, Names = NULL, Sep = ":") 
# {
#     # require(GenomicRanges)
#     if (is.null(Names)) {
#         Names <- paste(Chrom, as.integer(Start), as.integer(End), sep = Sep)
#     }
#     if (is.null(Strand)) {
#         Strand <- rep("*", length(Chrom))
#     }
#     Object <- GenomicRanges::GRanges(seqnames = Rle(Chrom), 
#         ranges = IRanges(Start, end = End, names = Names), 
#         strand = Rle(strand(Strand)))
#     return(Object)
# }
# Transcript_ranges <- MakeGRangesObject(Chrom = ExonicDataFrame[,"Chromosome"], Start = ExonicDataFrame[,"ExonStart"], 
# 	End = ExonicDataFrame[,"ExonEnd"], Strand =  ExonicDataFrame[,"Strand"], Names = paste(ExonicDataFrame[,"TxName"],ExonicDataFrame[,"ExonNumber"],sep="-"))
# elementMetadata(Transcript_ranges)[,"gene"] <- ExonicDataFrame[,"GeneName"]
# elementMetadata(Transcript_ranges)[,"transcript"] <- ExonicDataFrame[,"TxName"]
# elementMetadata(Transcript_ranges)[,"exon"] <- paste(ExonicDataFrame[,"TxName"],ExonicDataFrame[,"ExonNumber"],sep = "-")


gff <- makeTxDbFromGFF("/home/hdhiman/Data/work/hd/ref/picr/ncbi/GCF_003668045.1_CriGri-PICR_genomic.gff")
plot <- function(data, sample, chr, start, stop){
				data.gr <- lapply(data, function(x) {
					 gr <- makeGRangesFromDataFrame(x, keep.extra.columns=TRUE)})
	gtrack <- GenomeAxisTrack()
	atrack <- GeneRegionTrack(gff, 
				cex.group=1,
				name = "Genes", 
				showID = TRUE, 
				geneSymbol = TRUE, 
				background.title="#FFFFFF",
				collapseTranscripts = TRUE,
				col = "#000000", fill = "#F6EFEE",
				col.axis = "#000000",
				col.frame = "#000000",
				fontcolor = "#000000",
				fontcolor.title = "#000000")
	p1p1 <- DataTrack(range = data.gr[[3]], data = data.gr[[3]]$value, type = "histogram", col.histogram = NA, name = paste(sample," Primer 1 (P1)",sep=""), fill.histogram="#7DB8E8", background.title="#FFFFFF", col.title="#000000", col.axis="#696969")
	p2p1 <- DataTrack(range = data.gr[[4]], data = data.gr[[4]]$value, type = "histogram", col.histogram = NA, name = paste(sample," Primer 2 (P1)",sep=""), fill.histogram="#7DB8E8", background.title="#FFFFFF", col.title="#000000", col.axis="#696969")
	p1p10 <- DataTrack(range = data.gr[[1]], data = data.gr[[1]]$value, type = "histogram", col.histogram = NA, name = paste(sample," Primer 1 (P10)",sep=""), fill.histogram="#4682B4", background.title="#FFFFFF", col.title="#000000", col.axis="#696969")
	p2p10 <- DataTrack(range = data.gr[[2]], data = data.gr[[2]]$value, type = "histogram", col.histogram = NA, name = paste(sample," Primer 2 (P10)",sep=""), fill.histogram="#4682B4", background.title="#FFFFFF", col.title="#000000", col.axis="#696969")
pdf(paste(sample,".pdf", sep=""))
plotTracks(list(gtrack, atrack, p1p1, p2p1, p1p10, p2p10), 
	from = start, to = stop, chromosome = chr, 
	labelPos="below", main=chr, cex.main=0.8, 
	col.main="#000000", fontface.main=2, stacking = "dense",
	transcriptAnnotation = "gene")
dev.off()
}



plot(sample.G9, "G9", "NW_020822533.1", 6590000, 6620000)
plot(sample.5G10, "5G10.1", "NW_020822425.1", 6430000, 6730000)
plot(sample.5G10, "5G10.2", "NW_020822407.1", 9100000, 9800000)
plot(sample.5G10, "5G10.3", "NW_020822603.1", 1, 2400000)
plot(sample.1E3, "1E3", "NW_020822529.1", 13500000, 13900000)
plot(sample.1C11, "1C11", "NW_020822464.1", 4490000, 4530000)
plot(sample.6A6, "6A6_site1", "NW_020822506.1", 18290000, 18330000)
plot(sample.6A6, "6A6_site2", "NW_020822426.1", 2110000, 2150000)
plot(sample.6H1, "6H1_site1", "NW_020822506.1", 18290000, 18330000)
plot(sample.6H1, "6H1_site2", "NW_020822426.1", 2110000, 2150000)


