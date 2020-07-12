require("stringr")
require("data.table")
require("GenomicRanges")
require("Gviz")
require("genomation")
require("GenomicFeatures")
options(ucscChromosomeNames=FALSE) 

act<-read.table("src2/Tp5_11_segments_active.bed", header=F)
	colnames(act)<-c("chr", "start", "end", "state")
	act$col<-"#FFFF00"
enh<-read.table("src2/Tp5_11_segments_enh.bed", header=F)
	colnames(enh)<-c("chr", "start", "end", "state")
	enh$col<-"#FF0066"
prom<-read.table("src2/Tp5_11_segments_prom.bed", header=F)
	colnames(prom)<-c("chr", "start", "end", "state")
	prom$col<-"#006600"
qui<-read.table("src2/Tp5_11_segments_qui.bed", header=F)
	colnames(qui)<-c("chr", "start", "end", "state")
	qui$col<-"#D9D9D9"
rep<-read.table("src2/Tp5_11_segments_rep.bed", header=F)
	colnames(rep)<-c("chr", "start", "end", "state")
	rep$col<-"#002060"


gff <- makeTxDbFromGFF("/home/hdhiman/Data/work/hd/ref/picr/ncbi/GCF_003668045.1_CriGri-PICR_genomic.gff")
gtrack <- GenomeAxisTrack()
atrack <- GeneRegionTrack(gff, 
				cex.group=1,
				name = "Genes", 
				showID = TRUE, 
				geneSymbol = TRUE, 
				background.title="#FFFFFF",
				transcriptAnnotation = "GENEID",
				col = "#000000", fill = "#FFFF00",
				col.axis = "#000000",
				col.frame = "#000000",
				fontcolor = "#000000",
				fontcolor.title = "#000000")
var <- DataTrack(range = "src2/all_snp-2kb_mval_gt0.7.bedGraph", type="histogram",  
                 col.histogram = NA, name = "VariationPeak", fill.histogram="#F98484", 
                 background.title="#FFFFFF", col.title="#000000", col.axis="#696969")
# displayPars(var) <- list(alpha.title = 1, alpha = 0.5)
exp <- DataTrack(range = "src2/exp/bg/C3234A_rep1_mval.bedGraph", type="histogram",   
                 col.histogram = NA, name = "ExpressionPeak", fill.histogram="#1eb2a6", 
                 background.title="#FFFFFF", col.title="#000000", col.axis="#696969")

anno<-c(1:11)
names(anno)<-c("polycomb", "qui", "heterochrom", "trans", "weak_genic_enh", "weak_enh", "active_enh_1", "active_enh_2", "active_prom", "flanking_tss_upstream", "flanking_tss_downstream")
					 
states<-rbind(act,enh,prom,qui,rep)
	states.gr <- makeGRangesFromDataFrame(states, keep.extra.columns=TRUE)	
	states.gr$anno<-names(anno[states.gr$state])				 

plotRegion<-function(chr,a,b,c,d){
	states.track <- AnnotationTrack(states.gr[seqnames(states.gr)==chr], name = "States",
									col=NA,background.title="#FFFFFF", col.title="#000000", col.axis="#696969")				 
	feature(states.track)<-states.gr[seqnames(states.gr)==chr]$anno

	ht <- HighlightTrack(trackList = list(atrack, states.track, var, exp),
	                     start = a, width = b,
	                     chromosome = chr)
	plotTracks(list(gtrack, ht), from=c, to=d,
								chromosome = chr, labelPos="below", cex.main=0.8, col.main="#000000", 
								fontface.main=2, stacking = "dense", collapseTranscripts="longest", 
								transcriptAnnotation = "symbol", 
								groupAnnotation="anno", sizes=c(0.6,0.6,1,1,1),
								trans="#FFFF00", weak_enh="#FF0066", active_enh_1="#FF0066", 
								weak_genic_enh="#FF0066", active_enh_2="#FF0066", 
								flanking_tss_downstream="#006600", active_prom="#006600", 
								flanking_tss_upstream="#006600", qui="#FFFFFF", polycomb="#002060", 
								heterochrom="#002060")
}

pdf("site_stable.pdf", height=3.5)
plotRegion("NW_020822533.1", 6601451, 25, 6559451,6633451)
dev.off()

pdf("site_unstable.pdf", height=3.5)
plotRegion("NW_020822638.1", 1900568, 1, 1200568,2300568)
dev.off()
