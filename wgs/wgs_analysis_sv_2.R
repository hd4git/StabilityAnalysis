require("StructuralVariantAnnotation")
require("VariantAnnotation")
require("GenomicRanges")
require("rtracklayer")
require("UpSetR")
require("ComplexHeatmap")
require("reshape2")
require("ggplot2")
require("stringr")
require("graph")
require("igraph")
require("RColorBrewer")
require("circlize")
require("plyr")

setwd("/data/heena/work/hd/janssen/sv_analysis")

# Load gene annotation 
bed<-read.table("../../ref/picr/ncbi/GCF_003668045.1_CriGri-PICR_genic.bed", header=T)
genes<-makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE)

# Load chromosome association with scaffold names
chr.ass<-read.table("/data/heena/work/hd/ref/picr/ncbi/picr.chromosome_assignment.txt", header=F)
names(chr.ass)<-c("refseq", "picr", "chr")

# Create color table for samples with descriptive names
names<-c("C3234A","PE24","G9","5G10","1E3","1C11","6A6","6H1")
desc<-c("Host", "Process.Evolved", "LowCopy.Stable", "LowCopy.Unstable", "MediumCopy.Stable", "MediumCopy.Unstable", "HighCopy.Stable", "HighCopy.Unstable")
colors<-c("#71ACBC", "#3288BD", "#ABDDA4", "#66C2A5", "#FDAE61", "#F46D43","#D53E4F","#9E0142")
info<-data.frame(desc,colors)
rownames(info)<-names

legend<-data.frame(colors)
rownames(legend)<-names

cols<-as.character(info$colors)
names(cols)<-info$desc

############################
####### Load SV calls ######
############################
for(i in names) {
	assign(paste("inv.",i,sep=""), readVcf(paste("inversions/wgs_",i,"_inv_delly.vcf.gz", sep=""), "cgr"))
	assign(paste("dup.",i,sep=""), readVcf(paste("duplications/wgs_",i,"_dup_manta.vcf.gz", sep=""), "cgr"))
	assign(paste("ins.",i,sep=""), readVcf(paste("insertions/wgs_",i,"_ins_manta.vcf.gz", sep=""), "cgr"))
	assign(paste("del.manta.",i,sep=""), readVcf(paste("deletions/wgs_",i,"_del_manta.vcf.gz", sep=""), "cgr"))
	assign(paste("del.lumpy.",i,sep=""), readVcf(paste("deletions/wgs_",i,"_del_lumpy.vcf.gz", sep=""), "cgr"))
	assign(paste("trans.manta.",i,sep=""), readVcf(paste("translocations/wgs_",i,"_bnd_manta.vcf.gz", sep=""), "cgr"))
	assign(paste("trans.lumpy.",i,sep=""), readVcf(paste("translocations/wgs_",i,"_bnd_lumpy.vcf.gz", sep=""), "cgr"))
}


for(i in names) {
	assign(paste(i,"_cvg",sep=""), read.table(paste("cvg/",i,"_cvg/",i,"_cvg.txt", sep=""), header=FALSE))
}


############################################
############ Get SV granges ################
############################################

getGranges<-function( inv, del1, del2, dup, ins, trans1, trans2) {
	gr.inv <- breakpointRanges(inv)
	gr.ins <- breakpointRanges(ins)
	gr.dup <- breakpointRanges(dup)

	gr.del.m <- breakpointRanges(del1)
	gr.del.l <- breakpointRanges(del2)
	ol.del<-findBreakpointOverlaps(gr.del.l, gr.del.m,maxgap=100,sizemargin=0.25,restrictMarginToSizeMultiple=0.5)
	gr.del<-gr.del.m[subjectHits(ol.del),]

	gr.m1 <- breakpointRanges(trans1)
	gr.m2 <- breakendRanges(trans1)
	gr.m <- sort(c(gr.m1, gr.m2))
 	gr.l1 <- breakpointRanges(trans2)
  	gr.l2 <- breakendRanges(trans2)	
	gr.l <- sort(c(gr.l1, gr.l2))
	ol<-findBreakpointOverlaps(gr.l, gr.m,maxgap=1000,sizemargin=0.25,restrictMarginToSizeMultiple=0.5)
	gr.trans<-gr.m[subjectHits(ol),]

	granges.list<-list(inversions=gr.inv, insertions=gr.ins, deletions=gr.del, duplications=gr.dup, translocations=gr.trans)
	return(granges.list)
}


granges.C3234A<-suppressWarnings(getGranges(inv.C3234A, del.manta.C3234A, del.lumpy.C3234A, dup.C3234A, ins.C3234A, trans.manta.C3234A, trans.lumpy.C3234A))
granges.PE24<-suppressWarnings(getGranges(inv.PE24, del.manta.PE24, del.lumpy.PE24, dup.PE24, ins.PE24, trans.manta.PE24, trans.lumpy.PE24))
granges.G9<-suppressWarnings(getGranges(inv.G9, del.manta.G9, del.lumpy.G9, dup.G9, ins.G9, trans.manta.G9, trans.lumpy.G9))
granges.5G10<-suppressWarnings(getGranges(inv.5G10, del.manta.5G10, del.lumpy.5G10, dup.5G10, ins.5G10, trans.manta.5G10, trans.lumpy.5G10))
granges.1E3<-suppressWarnings(getGranges(inv.1E3, del.manta.1E3, del.lumpy.1E3, dup.1E3, ins.1E3, trans.manta.1E3, trans.lumpy.1E3))
granges.1C11<-suppressWarnings(getGranges(inv.1C11, del.manta.1C11, del.lumpy.1C11, dup.1C11, ins.1C11, trans.manta.1C11, trans.lumpy.1C11))
granges.6A6<-suppressWarnings(getGranges(inv.6A6, del.manta.6A6, del.lumpy.6A6, dup.6A6, ins.6A6, trans.manta.6A6, trans.lumpy.6A6))
granges.6H1<-suppressWarnings(getGranges(inv.6H1, del.manta.6H1, del.lumpy.6H1, dup.6H1, ins.6H1, trans.manta.6H1, trans.lumpy.6H1))

############ Remove repeated annotations to get correct counts of SVs
sv<-c("insertions", "inversions", "duplications", "deletions")
granges2.C3234A<-c(lapply(granges.C3234A[names(granges.C3234A) %in% sv], function(x) {x[grep("_bp1", names(x))]}), "translocations"=granges.C3234A[["translocations"]])
granges2.PE24<-c(lapply(granges.PE24[names(granges.PE24) %in% sv], function(x) {x[grep("_bp1", names(x))]}), "translocations"=granges.PE24[["translocations"]])
granges2.G9<-c(lapply(granges.G9[names(granges.G9) %in% sv], function(x) {x[grep("_bp1", names(x))]}), "translocations"=granges.G9[["translocations"]])
granges2.5G10<-c(lapply(granges.5G10[names(granges.5G10) %in% sv], function(x) {x[grep("_bp1", names(x))]}), "translocations"=granges.5G10[["translocations"]])
granges2.1E3<-c(lapply(granges.1E3[names(granges.1E3) %in% sv], function(x) {x[grep("_bp1", names(x))]}), "translocations"=granges.1E3[["translocations"]])
granges2.1C11<-c(lapply(granges.1C11[names(granges.1C11) %in% sv], function(x) {x[grep("_bp1", names(x))]}), "translocations"=granges.1C11[["translocations"]])
granges2.6A6<-c(lapply(granges.6A6[names(granges.6A6) %in% sv], function(x) {x[grep("_bp1", names(x))]}), "translocations"=granges.6A6[["translocations"]])
granges2.6H1<-c(lapply(granges.6H1[names(granges.6H1) %in% sv], function(x) {x[grep("_bp1", names(x))]}), "translocations"=granges.6H1[["translocations"]])


###################################
###### Get common SVs #############
###################################

common.sv<-function(x){
	base <- GenomicRanges::reduce(c(granges2.C3234A[[x]], 
									granges2.PE24[[x]], 
									granges2.G9[[x]], 
									granges2.5G10[[x]], 
									granges2.1E3[[x]],
									granges2.1C11[[x]],
									granges2.6A6[[x]],
									granges2.6H1[[x]]))
	z <- matrix(0,length(base),8)
	z[subjectHits(findOverlaps(granges.C3234A[[x]], base)),1] <- 1
	z[subjectHits(findOverlaps(granges.PE24[[x]], base)),2] <- 1
	z[subjectHits(findOverlaps(granges.G9[[x]], base)),3] <- 1
	z[subjectHits(findOverlaps(granges.5G10[[x]], base)),4] <- 1
	z[subjectHits(findOverlaps(granges.1E3[[x]], base)),5] <- 1
	z[subjectHits(findOverlaps(granges.1C11[[x]], base)),6] <- 1
	z[subjectHits(findOverlaps(granges.6A6[[x]], base)),7] <- 1
	z[subjectHits(findOverlaps(granges.6H1[[x]], base)),8] <- 1
	colnames(z) <- paste("gr.", 1:8)
	mcols(base) <- z
	base.df<-as.data.frame(base)
	base.df.same<-base.df[which(rowSums(base.df[,6:13])==8),]
	return(base.df.same)
}

type<-c("inversions", "insertions", "deletions", "duplications", "translocations")
for(i in type){
	assign(paste(i,"common", sep="."), common.sv(i))
}

common.sv.list<-list(makeGRangesFromDataFrame(inversions.common),
					makeGRangesFromDataFrame(insertions.common),
					makeGRangesFromDataFrame(deletions.common),
					makeGRangesFromDataFrame(duplications.common),
					makeGRangesFromDataFrame(translocations.common))
names(common.sv.list)<-type

##################################
###### Remove common #############
##################################

rmv.common<- function(x, y) {
	ol <- findOverlaps(x, y)
	Seq <- seq_along(x)
	return(x[!(Seq %in% queryHits(ol))])
}

master<-list(granges2.C3234A, granges2.PE24, granges2.G9, granges2.5G10, granges2.1E3, granges2.1C11, granges2.6A6, granges2.6H1)
names(master)<-names
names(type)<-type
master.f<-lapply(master, function(y) {lapply(type, function(x) {rmv.common(y[[x]], common.sv.list[[x]])})})

####################################################################
###### Assign chromosome annotation to create master file ##########
####################################################################

assign.chr2<-function(data){
		a<-chr.ass$chr[match(data$seqnames,chr.ass$refseq)]
		c<-data.frame(chr.ref=a)
		return(c)
	}

master.f2<-lapply(master.f, function(y) {lapply(type, function(x) {
	x.df<-as.data.frame(y[[x]], row.names=NULL)
	chr.info<-assign.chr2(x.df)
	cbind(x.df, chr.info)
	}
	)})



######################################
###### Get counts of SVs #############
######################################

counts<-do.call(rbind, lapply(master.f, function(sample){ vapply(sample, length, 1)}))
counts.m<-melt(counts)
colnames(counts.m)<-c("Sample", "SV.type", "Frequency")
counts.m$desc<-factor(info[counts.m$Sample,1], levels=as.vector(desc))

pdf(file="sv_counts_cor.pdf")
ggplot(data=counts.m, aes(x=desc, y=Frequency, fill=SV.type)) + 
	geom_bar(stat="identity", position="stack") +
	scale_fill_brewer(palette="Reds") + 
	labs(x="Cell lines", y="Frequency", fill="Structural variant") +
	theme_bw() +
	theme(legend.position="bottom") +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()

#####################################################
###### Print all SVs for all cell lines #############
#####################################################

samples.list.df<-lapply(master.f, function(x) {ldply(x, data.frame, .id=NULL)})
order<-c(1:3,11,9,5:8,10,4,12:17)

lapply(seq_along(samples.list.df), function(x) {
		temp<-samples.list.df[[x]][,order]
		temp[is.na(temp)]<-"."
		write.table(temp, paste("/var/www/html/jbrowse/data/sv_vcf/df", names(temp[x]),"sv.bed", sep="."),sep="\t",quote=F,row.names=F,col.names=F)
})

################################################
################## Translocations ##############
################################################

trans.list<-lapply(master.f, function(x) { x[["translocations"]]})

split.alt<-function(x){
	ref<-as.vector(seqnames(x))
	alt<-as.data.frame(str_split_fixed(x$ALT, "[\\[\\]\\:]", 4))[,2]
	trans.df<-data.frame(ref=ref, alt=alt)
	return(trans.df)
}


trans.df.list<-lapply(trans.list, function(x){split.alt(x)})

assign.chr<-function(x){
		a<-chr.ass$chr[match(x$ref,chr.ass$refseq)]
		b<-chr.ass$chr[match(x$alt,chr.ass$refseq)]
		c<-data.frame(chr.ref=a, chr.alt=b)
		return(c)
	}


for(i in names){
	info<-assign.chr(trans.df.list[[i]])
	trans.df.list[[i]]<-cbind(trans.df.list[[i]],info)
	trans.df.list[[i]]<-trans.df.list[[i]][which(trans.df.list[[i]]$chr.ref != "unplaced" & trans.df.list[[i]]$chr.alt != "unplaced"),]
	trans.df.list[[i]]$type=paste(trans.df.list[[i]]$chr.ref, trans.df.list[[i]]$chr.alt, sep=":")
	trans.df.list[[i]][which(trans.df.list[[i]]$chr.ref==trans.df.list[[i]]$chr.alt),]$type="intra"
}

intra<-as.data.frame(unlist(lapply(trans.df.list, function(x) {dim(x[which(x$type=="intra"), ])[1]})))
colnames(intra)<-"counts"
intra<-cbind(intra, desc)

intra.chr<-ldply(lapply(trans.df.list, function(x) {as.data.frame(table(x[which(x$type=="intra" ), 3]))}), data.frame)
colnames(intra.chr)<-c("Cell.lines", "chr", "Freq")
intra.chr$desc<-intra[intra.chr$Cell.lines,2]
intra.chr$color<-legend[intra.chr$Cell.lines,]
intra.chr<-intra.chr[which(intra.chr$chr != "unplaced"),]
intra.chr$desc<-factor(intra.chr$desc, levels=as.vector(desc))

size<-read.table("/home/hdhiman/Data/work/hd/ref/picr/ncbi/GCF_003668045.1_CriGri-PICR_genomic.fa.fai", header=F)
size.chr<-merge(size, chr.ass, by.x='V1', by.y='refseq')
size.list <- split(size.chr$V2, size.chr$chr)
size.each<-as.data.frame(vapply(size.list, sum, 1)) 
size.each$chr<-rownames(size.each)
colnames(size.each)<-c("size", "chr")
size.each<-size.each[which(size.each$chr!="unplaced"),]

intra.chr.size<-merge(intra.chr, size.each, by="chr")
intra.chr.size<-intra.chr.size[order(intra.chr.size$size),]
intra.chr.size$norm<-(intra.chr.size$Freq/intra.chr.size$size)*100000000

counts.master.f2<-lapply(master.f2, function(x){
	info<-sapply(x, function(y) {
	table(y$chr.ref)
	})
	info<-as.data.frame(info)
	sum<-rowSums(info)
	cbind(info, sum)	
})


get.norm<-function(x){
	x<-x[which(rownames(x)!="unplaced"),]
	x$size<-size.each$size[match(rownames(x), size.each$chr)]
	x$norm<-(x$sum/x$size)*1000000
	return(x$norm)
}

master.f.norm<-do.call(rbind,lapply(counts.master.f2, function(x){get.norm(x)}))
colnames(master.f.norm)<-size.each$chr
mean<-colMeans(master.f.norm)
sd<-colSds(master.f.norm)
master.dist<-data.frame(Mean=mean, SD=sd)
master.dist$chr<-rownames(master.dist)
master.norm.info<-rbind(master.f.norm, t(master.dist))

pdf(file="norm.counts.pdf", height=3, width=3.5)
	ggplot(data=master.dist, aes(x=factor(master.dist$chr, levels=c(1:10,"X")), y=Mean)) + 
	geom_point(aes(fill=Mean, size=SD), color="#000000", alpha=0.7, shape=21) +
	scale_fill_distiller(palette="YlOrRd", direction=1) +
	scale_size_continuous(range=c(1,5)) +
	labs(x="Chromosomes", y="Mean of frequencies\n normalized by size") +
	theme_bw() 
dev.off()
	# theme(legend.position="bottom")

#################################################
######## Intra-chromosomal translocations #######
#################################################

pdf(file="intra.chr.trans.corr.pdf")
ggplot(data=intra.chr, aes(x=factor(intra.chr$chr, levels=c(1:10,"X") ), y=Freq, fill=desc)) + 
geom_bar(stat="identity", position = "dodge") +
labs(x="Chromosomes")+
theme_minimal() +
scale_fill_manual(values=cols, limits=desc) + 
ylim(-100,175) +
theme( axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")) +
coord_polar(start = 0)
dev.off()


pdf(file="intra.chr.trans.corr2.pdf")
ggplot(data=intra.chr.size, aes(x=factor(intra.chr.size$chr, levels=c(1:10,"X") ), y=norm, fill=desc)) + 
geom_bar(stat="identity", position = "dodge") +
labs(x="Chromosomes")+
theme_minimal() +
scale_fill_manual(values=cols, limits=desc) + 
ylim(-100,175) +
theme( axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-2,4), "cm")) +
coord_polar(start = 0)
dev.off()

pdf(file="intra.chr.trans.norm.violin.pdf")
ggplot(data = intra.chr.size, aes(x=factor(intra.chr.size$chr, levels=c(1:10,"X","unplaced")), y=norm)) + 
	geom_violin() +
	geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3, aes(fill = desc), position=position_jitter(0.2)) +
	scale_fill_manual(values=cols) + 
	# coord_flip() +
	labs(x="Chromosomes", y="Frequency normalized by size", fill="Cell lines")+
	theme_bw()
dev.off()

#################################################
######## Inter-chromosomal translocations #######
#################################################

chr<-c(1:10,"X")
sel<-NULL
for(i in chr) {
	for(j in chr) {
		if(j>i) {
			n<-paste(i,j,sep=":") 
			sel<-c(sel,n)
		}
	}
}

# ggplot(as.data.frame(table(trans.df.list[[1]][which(trans.df.list[[1]]$type %in% sel),]$type)), aes(x=reorder(Var1, -Freq), y=Freq)) + geom_bar(stat="Identity") + theme_bw()
trans.freq.list<-lapply(trans.df.list, function(x) {
	freq<-as.data.frame(table(x[which(x$type %in% sel),]$type))
	freq2<-cbind(freq, as.data.frame(str_split_fixed(freq$Var1, ":", 2)))
	colnames(freq2)<-c("Alt", "weight", "from", "to")
	df<-freq2[,c(3,4,2)]
	# adj.m<-dcast(df, from~to, value.var="weight",fill=0)
	# row.names(adj.m)<- adj.m[,1]
	# adj.m<- as.matrix(adj.m[,-1])
	return(df)
})

grid.col=brewer.pal(n = 11, name = "Spectral")
pdf(file="interchromosomal_frequency.corr.pdf")
par(mfrow = c(2, 4))
chordDiagram(trans.freq.list[[1]], grid.col=grid.col)
title("Host")
chordDiagram(trans.freq.list[[2]], grid.col=grid.col)
title("Process Evolved")
chordDiagram(trans.freq.list[[3]], grid.col=grid.col)
title("Low copy, Stable")
chordDiagram(trans.freq.list[[4]], grid.col=grid.col)
title("Low copy, Untable")
chordDiagram(trans.freq.list[[5]], grid.col=grid.col)
title("Medium copy, Stable")
chordDiagram(trans.freq.list[[6]], grid.col=grid.col)
title("Medium copy, Untable")
chordDiagram(trans.freq.list[[7]], grid.col=grid.col)
title("High copy, Stable")
chordDiagram(trans.freq.list[[8]], grid.col=grid.col)
title("High copy, Untable")
dev.off()




################################################
############ Overlap with genes ################
################################################
alias<-str_split_fixed(genes$Alias, ",", 3)
genes$entrezID<-str_split_fixed(genes$Alias, ",", 3)[,1]
genes$type<-str_split_fixed(genes$Alias, ",", 3)[,3]

intersectGene<-function(inv, del, dup, ins, trans) {

	inv.gene.ol<-findOverlaps(inv, genes, type="any", maxgap=2000)
	gene.list.inv<-unique(genes[subjectHits(inv.gene.ol),]$entrezID)
	
	ins.gene.ol<-findOverlaps(ins, genes, type="any", maxgap=2000)
	gene.list.ins<-unique(genes[subjectHits(ins.gene.ol),]$entrezID)

	del.gene.ol<-findOverlaps(del, genes, type="any", maxgap=2000)
	gene.list.del<-unique(genes[subjectHits(del.gene.ol),]$entrezID)

	dup.gene.ol<-findOverlaps(dup, genes, type="any", maxgap=2000)
	gene.list.dup<-unique(genes[subjectHits(dup.gene.ol),]$entrezID)
	
	trans.gene.ol<-findOverlaps(trans, genes, type="any", maxgap=2000)
	gene.list.trans<-unique(genes[subjectHits(trans.gene.ol),]$entrezID)

	gene.ol.list<-list(inversions=gene.list.inv, insertions=gene.list.ins, deletions=gene.list.del, duplications=gene.list.dup, translocations=gene.list.trans)
	return(gene.ol.list)
}


gene.ol.list.G9<-suppressWarnings(intersectGene(master.f[["G9"]]$inversions, master.f[["G9"]]$insertions, master.f[["G9"]]$deletions, master.f[["G9"]]$duplications, master.f[["G9"]]$translocations))
gene.ol.list.5G10<-suppressWarnings(intersectGene(master.f[["5G10"]]$inversions, master.f[["5G10"]]$insertions, master.f[["5G10"]]$deletions, master.f[["5G10"]]$duplications, master.f[["5G10"]]$translocations))
gene.ol.list.1E3<-suppressWarnings(intersectGene(master.f[["1E3"]]$inversions, master.f[["1E3"]]$insertions, master.f[["1E3"]]$deletions, master.f[["1E3"]]$duplications, master.f[["1E3"]]$translocations))
gene.ol.list.1C11<-suppressWarnings(intersectGene(master.f[["1C11"]]$inversions, master.f[["1C11"]]$insertions, master.f[["1C11"]]$deletions, master.f[["1C11"]]$duplications, master.f[["1C11"]]$translocations))
gene.ol.list.6A6<-suppressWarnings(intersectGene(master.f[["6A6"]]$inversions, master.f[["6A6"]]$insertions, master.f[["6A6"]]$deletions, master.f[["6A6"]]$duplications, master.f[["6A6"]]$translocations))
gene.ol.list.6H1<-suppressWarnings(intersectGene(master.f[["6H1"]]$inversions, master.f[["6H1"]]$insertions, master.f[["6H1"]]$deletions, master.f[["6H1"]]$duplications, granges.6H1$translocations))
gene.ol.list.C3234A<-suppressWarnings(intersectGene(master.f[["C3234A"]]$inversions, master.f[["C3234A"]]$insertions, master.f[["C3234A"]]$deletions, master.f[["C3234A"]]$duplications, master.f[["C3234A"]]$translocations))
gene.ol.list.PE24<-suppressWarnings(intersectGene(master.f[["PE24"]]$inversions, master.f[["PE24"]]$insertions, master.f[["PE24"]]$deletions, master.f[["PE24"]]$duplications, master.f[["PE24"]]$translocations))

gene.ol.names<-list(unique(unlist(gene.ol.list.C3234A)), 
					unique(unlist(gene.ol.list.PE24)),
					unique(unlist(gene.ol.list.G9)), 
					unique(unlist(gene.ol.list.5G10)), 
					unique(unlist(gene.ol.list.1E3)), 
					unique(unlist(gene.ol.list.1C11)),
					unique(unlist(gene.ol.list.6A6)),
					unique(unlist(gene.ol.list.6H1)))
names(gene.ol.names)<-names

common<-Reduce(intersect,gene.ol.names)
length(common)

check<-function(x) {unique(unlist(x))[which(!unique(unlist(x)) %in% common)]}

gene.ol.count<-data.frame( 	lengths(lapply(gene.ol.list.C3234A, check)),
							lengths(lapply(gene.ol.list.PE24, check)),
							lengths(lapply(gene.ol.list.G9, check)), 
							lengths(lapply(gene.ol.list.5G10, check)), 
							lengths(lapply(gene.ol.list.1E3, check)), 
							lengths(lapply(gene.ol.list.1C11, check)),
							lengths(lapply(gene.ol.list.6A6, check)),
							lengths(lapply(gene.ol.list.6H1, check)))

colnames(gene.ol.count)<-names
unique.total<-lengths(lapply(gene.ol.names,check))
gene.ol.count<-rbind(gene.ol.count, unique.total)
rownames(gene.ol.count)<-c("inversions","insertions", "deletions", "duplications", "translocations", "unique.total")


gene.ol.count.m<-melt(t(as.matrix(gene.ol.count)))
colnames(gene.ol.count.m)<-c("Sample", "SV.type", "Frequency")
gene.ol.count.m$desc<-factor(intra[gene.ol.count.m$Sample,2], levels=as.vector(desc))

pdf(file="sv_gene_counts_woCommon_corr.pdf")
ggplot(data=gene.ol.count.m, aes(x=desc, y=Frequency, fill=SV.type)) + 
	geom_bar(stat="identity", position="dodge") +
	scale_fill_brewer(palette="Oranges") + 
	labs(x="Cell lines", y="Frequency", fill="Structural variant") +
	theme_bw() +
	theme(legend.position="bottom") +
	theme(axis.text.x = element_text(angle = 45, hjust=1)) 
dev.off()



gene.ol.names2<-lapply(gene.ol.names, check)
names(gene.ol.names2)<-desc
m<-make_comb_mat(gene.ol.names2)
# pdf(file="sv.genes.upsetPlot.pdf", onefile=T, paper='A4r')
# UpSet(m[comb_size(m) >= 5],
# 	set_order = desc, 
# 	pt_size = unit(2, "mm"), 
# 	lwd = 0.8,
# 	comb_col = brewer.pal(n=7, "Set2")[comb_degree(m)], 
# 	right_annotation = upset_right_annotation(m, gp = gpar(fill = colors)),
#     left_annotation = NULL,
#     show_row_names = FALSE
#    )
# dev.off()

pdf(file="sv.genes.upsetPlot.corr.pdf", onefile=T, paper='A4r')
UpSet(m,
	set_order = desc, 
	pt_size = unit(1, "mm"), 
	lwd = 0.8,
	comb_col = brewer.pal(n=7, "Set2")[comb_degree(m)], 
	left_annotation = rowAnnotation(
    "Set size" = anno_barplot(set_size(m), 
        border = FALSE, 
        gp = gpar(fill = colors), 
        width = unit(3, "cm")
    )), 
    right_annotation = NULL,
    show_row_names = FALSE
    )
dev.off()

combs<-make_comb_mat(gene.ol.names2)
comb_name(combs)

distinct.C3234A<-extract_comb(combs,"10000000")
distinct.PE24<-extract_comb(combs,"01000000")
distinct.G9<-extract_comb(combs,"00100000")
distinct.5G10<-extract_comb(combs,"00010000")
distinct.1E3<-extract_comb(combs,"00001000")
distinct.1C11<-extract_comb(combs,"00000100")
distinct.6A6<-extract_comb(combs,"00000010")
distinct.6H1<-extract_comb(combs,"00000001")

distinct.unstable<-c(extract_comb(combs,"00010101"),extract_comb(combs,"00010100"), extract_comb(combs,"00010001"), extract_comb(combs,"00000101"))
distinct.stable<-c(extract_comb(combs,"00101010"),extract_comb(combs,"00101000"), extract_comb(combs,"00100010"), extract_comb(combs,"00001010"))
distinct.unfavorable<-c(extract_comb(combs,"00010101"),extract_comb(combs,"00010111"),extract_comb(combs,"00010100"), extract_comb(combs,"00010001"), extract_comb(combs,"00000101"))

write.table(distinct.unstable,"distinct.unstable.txt",sep="\t",quote=F,row.names=F,col.names=F)
write.table(distinct.stable,"distinct.stable.txt",sep="\t",quote=F,row.names=F,col.names=F)




