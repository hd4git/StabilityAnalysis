library("DESeq2")
library("dplyr")
library("ggplot2")
library("vsn")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("genefilter")
library("apeglm")
library("ggpubr") 
library("UpSetR") 
library("GenomicRanges")

#Load the counts data
counts<-read.table("counts_gene_all.txt", header=T)

ID<-c("G9.rep1", "G9.rep2", "G9.rep3", "G10.rep1", "G10.rep2", "G10.rep3", "E3.rep1", "E3.rep2", "E3.rep3", "C11.rep1", "C11.rep2", "C11.rep3", "A6.rep1", "A6.rep2", "A6.rep3", "H1.rep1", "H1.rep2", "H1.rep3", "C1835A.rep1", "C1835A.rep2", "C1835A.rep3", "C3234A.rep1", "C3234A.rep2", "C3234A.rep3",  "PE24.rep1", "PE24.rep2", "PE24.rep3")
phase<-c(rep("G9",3),rep("G10",3),rep("E3",3),rep("C11",3),rep("A6",3),rep("H1",3),rep("C1835A",3),rep("C3234A",3),rep("PE24",3))
level<-c(rep("low",6),rep("med",6),rep("high",6),rep("host",9))
level2<-c(rep("lower",12),rep("higher",6),rep("host",9))
stability<-c(rep("stable",3),rep("unstable",3),rep("stable",3),rep("unstable",3),rep("stable",3),rep("unstable",3),rep("host",9))
type<-c(rep("exp",18), rep("host", 9))
type2<-c(rep("exp",18), rep("host2", 3),rep("host1", 6))
condition<-c(rep("exp",18),rep("host.C1835A",3),rep("host.C3234A",3),rep("host.PE24",3))
hosts<-c(rep("host1", 18), rep("host2", 3), rep("host1", 6))
favorability<-c(rep("fav",3),rep("unfav",3),rep("fav",3),rep("unfav",3),rep("unfav",3),rep("unfav",3),rep("host", 9))
conds <- as.factor(c(rep("LowCopy.Stable", 3), rep("LowCopy.Unstable", 3),
					rep("MediumCopy.Stable", 3), rep("MediumCopy.Unstable", 3),
					rep("HighCopy.Stable", 3), rep("HighCopy.Unstable", 3),
					rep("ATCC.Host", 3), rep("Horizon.Host", 3),rep("Process.Evolved", 3)))

sample_ann <- data.frame(ID, phase,level,level2, stability,type, type2, condition, hosts, favorability, conds)

colnames(counts)<-paste(conds,rep(c("1","2","3"), 9),sep=".")

################################################################  Part 1  ################################################################################


#######################################
######## Overview of all samples ######
#######################################

#Create input for DE analysis
batch <- DESeqDataSetFromMatrix( countData = counts,colData  = sample_ann, design = ~ phase)
batch <- estimateSizeFactors(batch)
batch <- estimateDispersions(batch)
#Filter lowly expressed genes
keep <- rowSums(counts(batch) >= 10) >= 3
batch.filter<- batch[keep,]
batch<-batch.filter

### MeanSD plot to see variance
ntd <- normTransform(batch)
vsd <- vst(batch)

## No normalization
pdf("out/meanSdPlot.pdf")
	meanSdPlot(assay(batch), ranks = FALSE)
dev.off()
## log Transform
pdf("out/meanSdPlot-ntd.pdf")
	meanSdPlot(assay(ntd), ranks = FALSE)
dev.off()
## VST normalization
pdf("out/meanSdPlot-vsd.pdf")
	meanSdPlot(assay(vsd), ranks = FALSE)
dev.off()


#Estimate library size factors and dispersion from the data and fit model to estimate Differential Expression 
DDS_batch <- DESeq(batch)

norm.countsbatch <- counts(DDS_batch, normalized=TRUE)
# conds = as.factor(sample_ann$phase)
cond_colours = c(rep("#ABDDA4",3),rep("#66C2A5",3), rep("#FDAE61",3), rep("#F46D43",3),  rep("#D53E4F",3), rep("#9E0142", 3), rep("#62D1EF",3),rep("#71ACBC",3), rep("#3288BD",3))
names(cond_colours)=conds

#### Plot PCA for top 500 variant genes
pcaData <- plotPCA(vsd, intgroup = "conds", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# pdf("out/pca-vsd_top500.pdf", width=9)
a<-ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(conds, level=unique(conds)))) +
  geom_point(size=5) +  
  scale_colour_manual(values=cond_colours) + 
  geom_text(label=rep(c("1","2","3"),9), color="#000000") +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"), 
  		y=paste0("PC2: ", percentVar[2], "% variance"),
  		color="Cell Lines") +
  theme_bw() +
  	theme(legend.position="bottom") 
# dev.off()

pcaData <- plotPCA(vsd[,1:18], intgroup = "conds", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# pdf("out/pca-vsd_top500_exp.pdf")
b<-ggplot(pcaData, aes(x = PC1, y = PC2, color = factor(conds[1:18], level=unique(conds[1:18])))) +
  geom_point(size=5) +  
  scale_colour_manual(values=cond_colours[1:18]) + 
  geom_text(label=rep(c("1","2","3"),6), color="#000000") +
  labs(x=paste0("PC1: ", percentVar[1], "% variance"), 
  		y=paste0("PC2: ", percentVar[2], "% variance"),
  		color="Cell Lines") +
  theme_bw() +
  	theme(legend.position="bottom") 
# dev.off()

### Plot heatmap with VST 
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$ID )
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
# ggsave("out/heatmap-vsd.pdf",device = "pdf",

rownames(sampleDistMatrix)<-paste(conds,rep(c("1","2","3"), 9),sep=".")
c<- as.grob(pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors))

figure <- ggarrange(a, c, 
                    labels = c("A", "B"),
                    ncol = 2, nrow = 1)

pdf("out/pca_heatmap.pdf", width=15)
figure
dev.off()

#### Plot PCA for all data

pca<-prcomp(assay(vsd), center=TRUE, scale=TRUE)
pdf(file="out/pca-vsd_all.pdf")
plot(pca$rotation, col=cond_colours, pch=19, cex=2, main= "PCA based on expression profiles with VST normalization")
text(pca$rotation, rep(c("1","2","3"),9), pos=3,cex=0.6)
legend("bottomright", col=unique(cond_colours), legend = unique(conds),pch = 20, bty='n', cex=.75)
dev.off()

pca<-prcomp(assay(vsd)[,1:18], center=TRUE, scale=TRUE)
pdf(file="out/pca-vsd_all_exp.pdf")
plot(pca$rotation, col=cond_colours, pch=19, cex=2, main= "PCA based on expression profiles with VST normalization")
text(pca$rotation, rep(c("1","2","3"),9), pos=3,cex=0.6)
legend("bottomright", col=unique(cond_colours), legend = unique(conds),pch = 20, bty='n', cex=.75)
dev.off()




# #Filter lowly expressed genes
# keep <- rowSums(counts(batch) >= 10) >= 3
# batch.filter<- batch[keep,]

fpm=fpm(batch, robust=TRUE)
keep2 <- rowSums(fpm >= 10) >= 3
fpm.filter<- fpm[keep,]
write.table(fpm.filter,"out/fpm.txt",sep="\t",quote=F,row.names=T,col.names=T)


################################################################  Part 2  ################################################################################

##############################
######## DE host vs exp ###
##############################
dds<-batch.filter
dds$type <- factor(dds$type, levels = c("host", "exp"))
design(dds) <- formula(~ type)

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

############################################################
### Using coef with lfcShrink for getting DEGs 
de_results<-lfcShrink(DDS_batch, coef="type_exp_vs_host", type="apeglm")
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/host.exp/host.exp.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
colnames(mat)<-conds
anno <- as.data.frame(colData(vsd)[, c("conds","condition")])
rownames(anno)<-paste(conds,rep(c("1","2","3"), 9),sep=".")

# conds1 <- unique(as.factor(sample_ann$phase))
# cond_colours1 <- c("#F2CCC3","#E78F8E", "#bde1dd", "#49bdb6",  "#B2DF8A", "#33A02C", "#eaadbd","#b88a9f", "#876880")
# names(cond_colours1)=conds1
conds2 <- unique(as.factor(sample_ann$condition))
cond_colours2 <- c("#b3b2b1", "#3288BD","#3288BD", "#3288BD" )
names(cond_colours2)=conds2
# color_list<-list(phase=cond_colours1, condition=cond_colours2)


# pdf(file="out/host.exp/heatmap_exp_host.pdf")
# 	pheatmap(mat, annotation_col = anno, 
# 					fontsize_row=10,  
# 					treeheight_row = 0, 
# 					treeheight_col = 0,
# 					cluster_rows=F, 
# 					cluster_cols=F)
# dev.off()

pdf(file="out/host.exp/heatmap_exp_host.pdf")
pheatmap(mat, 
	annotation_col = anno,
	annotation_colors = list(conds=cond_colours, condition=cond_colours2), 
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlGn")))(100),
 	fontsize_row=7,  
	treeheight_row = 0, 
	treeheight_col = 0,
	show_rownames = F,
	cluster_rows=F, 
	cluster_cols=F, 
	angle_col=45)
dev.off()

ggsave("out/host.exp/MAplot-DE_host_exp.pdf", ggmaplot(de_results, main = expression("Host Vs Expressing"),
   fdr = 0.01, fc = 1.5, size = 0.4,
   palette = c("#B31B21", "#1465AC", "darkgray"),
   genenames = as.vector(de_results$name),
   legend = "top", top = 60,
   font.label = 6,
   font.legend = "bold",
   font.main = "bold",
   ggtheme = ggplot2::theme_minimal()))


############################################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("type","exp","host"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/host.exp/host_exp.rnk",sep="\t",quote=F,row.names=T,col.names=F)


#####################################
######## DE host(w/o ATCC) vs exp ###
#####################################

dds<-batch.filter
dds$type2 <- factor(dds$type2, levels = c("host1","host2", "exp"))
design(dds) <- formula(~ type2)

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

############################################################
### Using coef with lfcShrink for getting DEGs 
de_results<-lfcShrink(DDS_batch, coef="type2_exp_vs_host1", type="apeglm")
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/host.exp/host1.exp.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/host.exp/heatmap_exp_host1.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()

############################################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("type2","exp","host1"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/host.exp/host1_exp.rnk",sep="\t",quote=F,row.names=T,col.names=F)

##########################################
######## DE host ATCC vs C3234A and PE ###
##########################################

############################################################
### Using coef with lfcShrink for getting DEGs 
de_results<-lfcShrink(DDS_batch, coef="type2_host2_vs_host1", type="apeglm")
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/host.exp/host1.host2.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/host.exp/heatmap_host2_host1.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()

############################################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("type2","host2","host1"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/host.exp/host1_host2.rnk",sep="\t",quote=F,row.names=T,col.names=F)


#####################################
######## DE host C3234A vs PE24 ###
#####################################

dds<-batch.filter
dds$condition <- factor(dds$condition, levels = c("host.C3234A","host.PE24", "host.C1835A", "exp"))
design(dds) <- formula(~ condition)

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

############################################################
### Using coef with lfcShrink for getting DEGs 
de_results<-lfcShrink(DDS_batch, coef="condition_host.PE24_vs_host.C3234A", type="apeglm")
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/host.exp/C3234A.PE24.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/host.exp/heatmap_C3234A_PE24.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()

############################################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("condition","host.PE24","host.C3234A"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/host.exp/C3234A_PE24.rnk",sep="\t",quote=F,row.names=T,col.names=F)


################################################################  Part 3  ################################################################################


#####################################
######## DE stable vs unstable ###
#####################################

dds<-batch.filter
dds$stability <- factor(dds$stability, levels = c("stable","unstable", "host"))
design(dds) <- ~ stability

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

############################################################
### Using coef with lfcShrink for getting DEGs 
de_results<-lfcShrink(DDS_batch, coef="stability_unstable_vs_stable", type="apeglm")
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/stable.unstable/stable.unstable.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/stable.unstable/heatmap_stable_unstable.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()

############################################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("stability","unstable","stable"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/stable.unstable/stable_unstable.rnk",sep="\t",quote=F,row.names=T,col.names=F)


#####################################
######## DE low vs med ###
#####################################
dds<-batch.filter
dds$level <- factor(dds$level, levels = c("low","med", "high", "host"))
design(dds) <- ~ level

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

############################################################
### Using coef with lfcShrink for getting DEGs 
de_results<-lfcShrink(DDS_batch, coef="level_med_vs_low", type="apeglm")
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/low.med/low.med.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/low.med/heatmap_low_med.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()

############################################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("level","med","low"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/low.med/low_med.rnk",sep="\t",quote=F,row.names=T,col.names=F)

#####################################
######## DE low vs high ###
#####################################
dds<-batch.filter
dds$level <- factor(dds$level, levels = c("low","med", "high", "host"))
design(dds) <- ~ level

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

############################################################
### Using coef with lfcShrink for getting DEGs 
de_results<-lfcShrink(DDS_batch, coef="level_high_vs_low", type="apeglm")
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/low.high/low.high.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/low.high/heatmap_low_high.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()

############################################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("level","high","low"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/low.high/low_high.rnk",sep="\t",quote=F,row.names=T,col.names=F)

#####################################
######## DE med vs high ###
#####################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("level","high","med"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/med.high/med_high.rnk",sep="\t",quote=F,row.names=T,col.names=F)



#####################################
######## DE low + med vs high ###
#####################################
dds<-batch.filter
dds$level2 <- factor(dds$level2, levels = c("lower","higher", "host"))
design(dds) <- ~ level2

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

############################################################
### Using coef with lfcShrink for getting DEGs 
de_results<-lfcShrink(DDS_batch, coef="level2_higher_vs_lower", type="apeglm")
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/low.med/higher.lower.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/low.med/heatmap_higher_lower.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()

############################################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("level2","higher","lower"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/low.med/higher_lower.rnk",sep="\t",quote=F,row.names=T,col.names=F)

################################################################  Part 4  ################################################################################

#Create input for DE analysis
dds<-batch.filter
dds$group <- factor(paste0(dds$level, dds$stability))
design(dds) <- ~ group

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

#####################################
######## DE G9 vs G10 ###
#####################################
de_results<-results(DDS_batch, contrast = c("group", "lowunstable", "lowstable"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585& padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/G9.G10/G9.G10.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/G9.G10/G9_G10.rnk",sep="\t",quote=F,row.names=T,col.names=F)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/G9.G10/heatmap_G9_G10.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()


#####################################
######## DE E3 vs C11 ###############
#####################################
de_results<-results(DDS_batch, contrast = c("group", "medunstable", "medstable"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585& padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/E3.C11/E3.C11.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/E3.C11/E3_C11.rnk",sep="\t",quote=F,row.names=T,col.names=F)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/E3.C11/heatmap_E3_C11.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()

#####################################
######## DE 6A6 vs 6H1 ###############
#####################################
de_results<-results(DDS_batch, contrast = c("group", "highunstable", "highstable"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585& padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/6A6.6H1/6A6.6H1.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/6A6.6H1/6A6_6H1.rnk",sep="\t",quote=F,row.names=T,col.names=F)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

pdf(file="out/6A6.6H1/heatmap_6A6_6H1.pdf")
pheatmap(mat, annotation_col = anno, fontsize_row=10,  treeheight_row = 0, treeheight_col = 0,cluster_rows=F, cluster_cols=F)
dev.off()


#####################################
######## DE fav vs unfav ###
#####################################

dds<-batch.filter
dds$favorability <- factor(dds$favorability, levels = c("fav","unfav", "host"))
design(dds) <- ~ favorability

DDS_batch <- DESeq(dds)
resultsNames(DDS_batch)

############################################################
### Using coef with lfcShrink for getting DEGs 
de_results<-lfcShrink(DDS_batch, coef="favorability_unfav_vs_fav", type="apeglm")
summary(de_results)

#Export sig DE list
sig_de_results <-subset(de_results,  abs(log2FoldChange)>0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)
write.table(as.data.frame(sig_de_results),file ="out/stable.unstable/fav.unfav.BH.0.01.FC.txt",sep="\t",quote=F,row.names=T,col.names=T)

topVarGenes <- rownames(sig_de_results)
nCounts <- counts(DDS_batch , normalized=TRUE)
mat<-as.matrix(log2(nCounts[topVarGenes, ]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])

ecm.genes <- read.table("fav_ecm.txt")
nCounts <- counts(DDS_batch , normalized=TRUE)
rownames(nCounts) <- toupper(rownames(nCounts))
mat<-as.matrix(log2(nCounts[which(rownames(nCounts) %in% ecm.genes$V1), 1:18]+1))
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("phase","condition")])


conds1 <- unique(as.factor(sample_ann$phase))
cond_colours1 <- c("#F2CCC3","#E78F8E", "#bde1dd", "#49bdb6",  "#B2DF8A", "#33A02C", "#eaadbd","#b88a9f", "#876880")
# brewer.pal(length(unique(conds1)),"Paired")[as.factor(conds1)]
names(cond_colours1)=conds1
conds2 <- unique(as.factor(sample_ann$condition))
cond_colours2 <- c("#b3b2b1", "#eaadbd","#b88a9f", "#876880" )
# brewer.pal(length(unique(conds2)),"Pastel2")[as.factor(conds2)]
names(cond_colours2)=conds2
color_list<-list(phase=cond_colours1, condition=cond_colours2)

# pdf(file="out/stable.unstable/heatmap_fav_unfav.pdf")
pheatmap(mat, 
	annotation_col = anno,
	annotation_colors = list(phase=cond_colours1, condition=cond_colours2), 
	color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlGn")))(100),
 	fontsize_row=7,  
	treeheight_row = 0, 
	treeheight_col = 0,
	cluster_rows=F, 
	cluster_cols=F, 
	angle_col=45)
dev.off()


#Plot Sig DE genes in MAplot
ggsave("out/stable.unstable/MAplot-DE_fav_unfav.pdf", ggmaplot(de_results, main = expression("Favorable Vs Unfavorable"),
   fdr = 0.01, fc = 1.5, size = 0.4,
   palette = c("#B31B21", "#1465AC", "darkgray"),
   genenames = as.vector(de_results$name),
   legend = "top", top = 60,
   font.label = 6,
   font.legend = "bold",
   font.main = "bold",
   ggtheme = ggplot2::theme_minimal()))

############################################################
### Using contrast for getting ranks for GSEA 
de_results<-results(DDS_batch, contrast = c("favorability","unfav","fav"),lfcThreshold=0, independentFiltering =T, pAdjustMethod="BH", alpha = 0.01)
#Export sig DE list
summary(de_results)

sig_de_results <-subset(de_results,  abs(log2FoldChange)> 0.585 & padj < 0.01)
sig_de_results <- sig_de_results[order(sig_de_results$log2FoldChange, decreasing=T),]
summary(sig_de_results)

res<-as.data.frame(de_results)
res<-res[order(res$stat,decreasing=TRUE),]
res2<-subset(res, select="stat")
rownames(res2)<-toupper(rownames(res2))
write.table(res2,file ="out/stable.unstable/fav_unfav.rnk",sep="\t",quote=F,row.names=T,col.names=F)

############################################################
# Bar plots for up and down regulated genes in fovorable/unfavorable comparisons

gsea.pos<-read.table("GSEA/fav_unfav_gsea_pos.txt", header=TRUE, sep="\t")
gsea.pos<-gsea.pos[order(gsea.pos$FDR.q.val),]
gsea.neg<-read.table("GSEA/fav_unfav_gsea_neg.txt", header=TRUE, sep="\t")
gsea.neg<-gsea.neg[order(gsea.neg$FDR.q.val),]

# gsea.pos<-gsea.pos[!grepl('BIOCARTA|PID', gsea.pos$NAME),]
gsea.pos<-gsea.pos[which(gsea.pos$FDR.q.val<=0.1),]

gsea.neg<-gsea.neg[which(gsea.neg$FDR.q.val<=0.1),]

gsea.neg$NAME <- factor(gsea.neg$NAME, levels = gsea.neg$NAME[order(gsea.neg$FDR.q.val, decreasing=TRUE)])
gsea.pos$NAME <- factor(gsea.pos$NAME, levels = gsea.pos$NAME[order(gsea.pos$FDR.q.val, decreasing=TRUE)])

pdf(file="out/stable.unstable/gsea_fav_unfav_neg.pdf", width=13,height=10)
ggplot(data=gsea.neg, aes(x=NAME, y=FDR.q.val, fill=NOM.p.val))+ 
geom_bar(stat="identity")+coord_flip() +
scale_fill_gradient(low="#A6D96A",high="#1A9641") + 
theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

pdf(file="out/stable.unstable/gsea_fav_unfav_pos.pdf", width=13,height=10)
gsea.pos$NAME <- factor(gsea.pos$NAME, levels = gsea.pos$NAME[order(gsea.pos$FDR.q.val, decreasing=TRUE)])
ggplot(data=gsea.pos, aes(x=NAME, y=FDR.q.val, fill=NOM.p.val))+ 
geom_bar(stat="identity")+coord_flip() +
scale_fill_gradient(low="#FDAE61",high="#D7191C") + 
theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

############################################################
# Bar plots for chromosome assignment to landing pads

fav<-read.table("favorable_chr.txt", header=F)
unfav<-read.table("unfavorable_chr.txt", header=F)
fav<-fav[match(vec, fav$V1),]
unfav<-unfav[match(vec, unfav$V1),]
pdf(file="favorable_chr.pdf")
ggplot(data=fav, aes(x=V1, y=V2)) + geom_bar(stat="identity", fill="#A6D96A") + labs(x="Chromosome", y="Frequency") + theme_minimal() 
dev.off()
pdf(file="unfavorable_chr.pdf")
ggplot(data=unfav, aes(x=V1, y=V2)) + geom_bar(stat="identity", fill="#FDAE61") + labs(x="Chromosome", y="Frequency") + theme_minimal() 
dev.off()

