require("factoextra")
require("ggfortify")
require("gridExtra")
require("ggrepel")
require("ggplot2")
require("ggpubr")

data<-read.table("integration_site_profile_new2.txt", header=TRUE, sep="\t")
# data<-transform(data, Expression.Peak.Distance= pmin(High.expression.bin.to.the.left, High.expression.bin.to.the.right))
cols<-c("Stable"="#71ACBC", 
		"Unstable"="#3288BD", 
		"LowCopy.Stable"="#ABDDA4", 
		"LowCopy.Stable.External"="#7DC673", 
		"LowCopy.Unstable"="#66C2A5",
		"MediumCopy.Stable"="#FDAE61", 
		"MediumCopy.Unstable"="#F46D43",
		"HighCopy.Stable"="#D53E4F",
		"HighCopy.Unstable"="#9E0142")

cols<-c("Stable"="#a6cee3", 
		"Unstable"="#1f78b4", 
		"LowCopy.Stable"="#b2df8a", 
		"LowCopy.Stable.External"="#cab2d6", 
		"LowCopy.Unstable"="#33a02c",
		"MediumCopy.Stable"="#fdbf6f", 
		"MediumCopy.Unstable"="#ff7f00",
		"HighCopy.Stable"="#fb9a99",
		"HighCopy.Unstable"="#e31a1c")

anno<-c("Stable", "Unstable", "LowCopy.Stable", "LowCopy.Stable.External", "LowCopy.Unstable", "MediumCopy.Stable", "MediumCopy.Unstable", "HighCopy.Stable", "HighCopy.Unstable")
data$desc<-factor(data$desc, levels=anno)

v<-ggplot(data = data, aes(x=factor(Stability), y=Variant.Profile)) + 
	geom_violin(trim = FALSE) +
	geom_boxplot(width=.1, color="grey")+ 
	geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, aes(fill = desc), position=position_jitter(w=0.1, h=0)) +
	scale_fill_manual(values=cols) + 
	# coord_flip() +
	labs(x="Stability", y="Variant profile at the site", fill="Cell lines")+
	theme_light()
e<-ggplot(data = data, aes(x=factor(Stability), y=Expression.Peak.Distance)) + 
	geom_violin(trim = FALSE) +
	geom_boxplot(width=.1, color="grey")+ 
	geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, aes(fill = desc), position=position_jitter(w=0.1, h=0)) +
	scale_fill_manual(values=cols) + 
	# coord_flip() +
	labs(x="Stability", y="Expression profile at the site", fill="Cell lines")+
	theme_light()
s<-ggplot(data = data, aes(x=factor(Stability), y=State)) + 
	geom_violin(trim = FALSE) +
	geom_boxplot(width=.1, color="grey")+ 
	geom_dotplot(binaxis='y', stackdir='center', dotsize=1, aes(fill = desc), position=position_jitter(w=0.1, h=0)) +
	scale_fill_manual(values=cols) + 
	# coord_flip() +
	labs(x="Stability", y="Chromatin state at the site", fill="Cell lines")+
	theme_light()
figure <- ggarrange(v, e, s, 
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1,   common.legend = TRUE, legend = "bottom")
pdf(file = "integration_site_profile_new.pdf", height=5, width=9)
figure
dev.off()








# v<-ggplot(data, aes(x=factor(Stability), y=Variant.Profile)) + 
# 	geom_violin(trim = FALSE)+ 
# 	geom_boxplot(width=.1, color="grey")+ 
# 	geom_point(aes(color = factor(desc)),  position=position_jitter(w=0.1, h=0), show.legend=T) + 
# 	scale_fill_manual(values=data$cols, limits=data$desc) + 
# 	xlab("Stability") + 
# 	ylab("Variant profile at the site") +
# 	theme_light() 
# 	dev.off()
# e<-ggplot(data, aes(x=factor(Stability), y=Expression.Peak.Distance)) + 
# 	geom_violin(trim = FALSE)+ 
# 	xlab("Stability") + 
# 	geom_boxplot(width=.1, color="grey")+ 
# 	ylab("Expression profile at the site") +
# 	geom_point(aes(color = factor(Favorability)),  position=position_jitter(w=0.1, h=0), show.legend=F) + 
# 	theme_light() 
# s<-ggplot(data, aes(x=factor(Stability), y=State)) + 
# 	geom_violin(trim = FALSE)+ xlab("Stability") + 
# 	geom_boxplot(width=.1, color="grey")+ 
# 	ylab("Chromatin state the site") +
# 	geom_point(aes(color = factor(Favorability)),  position=position_jitter(w=0.1, h=0)) + 
# 	theme_light() 
# # g<-ggplot(data, aes(x=factor(Stability), y=In.gene)) + geom_violin(trim = FALSE)+ xlab("Stability") + geom_boxplot(width=.1, color="grey")+ ylab("Overlap with a genic loci") +geom_point(aes(color = factor(Favorability)),  position=position_jitter(w=0.1,h=0.05)) + theme_light()
# figure <- ggarrange(v, e, s, 
#                     labels = c("A", "B", "C"),
#                     ncol = 3, nrow = 1,   common.legend = TRUE, legend = "bottom")
# # png(filename = "integration_site_profile_new.png")
# # figure
# # dev.off()
# pdf(file = "integration_site_profile_new.pdf", height=4)
# figure
# dev.off()













# data<-read.table("complete_analysis.txt", header=TRUE, sep="\t")
# df<-data[,c(6:19)]
# rownames(df)<-make.names(data$stability, unique=TRUE)
# fviz_nbclust(df, kmeans, method = "wss")+ labs(subtitle = "Elbow method")### Find optimal number of clusters by Elbow method. 
# autoplot(kmeans(df, 4), data = df, label=TRUE)

data<-read.table("complete_analysis2.txt", header=TRUE, sep="\t")

data["phenotype"]<-rep('1', times=44)
data$phenotype[data$Favorability.numeric!=3]<-0
data$phenotype<-as.numeric(data$phenotype)

data["stab"]<-rep('1', times=44)
data$stab[data$Favorability.numeric<0]<-0
data$stab<-as.numeric(data$stab)
data<-transform(data, min.expression.dis= pmin(High.expression.bin.to.the.left, High.expression.bin.to.the.right))
data<-transform(data, min.variant.dis= pmin(High.variant.bin.to.the.left, High.variant.bin.to.the.right))

p1<-ggplot(data, aes(x=In.gene, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + display.brewer.pal(9, "Set1") + ggtitle("Presence within a gene") + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() 
p2<-ggplot(data, aes(x=Variant.Profile, y=Favorability.numeric, color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + display.brewer.pal(9, "Set1")+ ggtitle("Variant profile at the site") + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15))
p3<-ggplot(data, aes(x=Expression.Profile, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + display.brewer.pal(9, "Set1")+ ggtitle("Expression profile at the site") + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15))
p4<-ggplot(data, aes(x=Chromatin.state, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + display.brewer.pal(9, "Set1")+ ggtitle("Chromatin state at the site")  + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15))
p5<-ggplot(data, aes(x=Chromatin.state.left, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + display.brewer.pal(9, "Set1")+ ggtitle("Chromatin state (left)")  + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15))
p6<-ggplot(data, aes(x=Chromatin.state.right, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F)+ display.brewer.pal(9, "Set1")+ ggtitle("Chromatin state (right)")  + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15))
p7<-ggplot(data, aes(x=state, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F)+ display.brewer.pal(9, "Set1")+ ggtitle("Chromatin state at/around the site")  + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15))
p8<-ggplot(data[1:24,], aes(x=Uniqueness.within.200bp, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F)+ display.brewer.pal(9, "Set1")+ ggtitle("Uniqueness within integration site")  + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15)) +theme(legend.position="bottom")
p9<-ggplot(data, aes(x=min.variant.dis, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + display.brewer.pal(9, "Set1")+ ggtitle("Distance from high variant bin")  + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15))
p10<-ggplot(data, aes(x=min.expression.dis, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F)+ display.brewer.pal(9, "Set1")+ ggtitle("Distance from high expression bin")  + theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15))
p11<-ggplot(data, aes(x=In.gene, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15)) +ggtitle("Presence within a gene") + display.brewer.pal(9, "Set1")+ theme(plot.title = element_text(size=6)) + theme(text = element_text(size=5)) + theme_light() + scale_shape_manual(values=seq(0,15)) +theme(legend.position="bottom")

pdf(file = "univar_relation_at_site_tg.pdf")
grid.arrange(p2, p3, p4,p5,p6, p7, nrow = 3)
grid.arrange( p8,p1, p9, p10,p11, nrow = 3, ncol=2)
dev.off()


data.f<-read.table("janssen_complete_analysis.txt", header=TRUE, sep="\t")
pdf(file = "univar_relation_fusion.pdf")
p1<-ggplot(data.f, aes(x=Structural.Variation , y=Favorability.numeric, shape=Favorability, color=Favorability)) + geom_point(size=3, position=position_jitter(h=0.15,w=0.15)) + ggtitle("Structural Variation") + theme(plot.title = element_text(size=8)) + theme(text = element_text(size=6)) + theme_light() + scale_shape_manual(values=seq(0,15))
p2<-ggplot(data.f, aes(x=TG.TG.Fusions, y=Favorability.numeric, shape=Favorability, color=Favorability)) + geom_point(size=3, position=position_jitter(h=0.15,w=0.15), show.legend=F) + ggtitle("Number of TG-TG fusions") + theme(plot.title = element_text(size=8)) + theme(text = element_text(size=6)) + theme_light() + scale_shape_manual(values=seq(0,15))
p3<-ggplot(data.f, aes(x=Percentage.Proper.orientation, y=Favorability.numeric, shape=Favorability, color=Favorability)) + geom_point(size=3, position=position_jitter(h=0.15,w=0.15), show.legend=F) + ggtitle("Percentage of proper orientation") + theme(plot.title = element_text(size=8)) + theme(text = element_text(size=6)) + theme_light() + scale_shape_manual(values=seq(0,15))
grid.arrange(p1, p2, p3, nrow = 3, ncol=2)
dev.off()

data2<-data[,-c(10:18,20,21)]
data.mat <- c("Field","good.mean","bad.mean","p.val.mean","good.median","bad.median","p.val.median")
for (x in 7:length(colnames(data2))){
	m1<-mean(as.matrix(data2[which(data2$stability=="stable" & data2$copy=="low"),x]))
	m2<-mean(as.matrix(data2[-(which(data2$stability=="stable" & data2$copy=="low")),x]))
	m3<-wilcox.test(max(m1,m2), min(m1,m2), paired = TRUE, alternative = "greater")$p.value

	m4<-median(as.matrix(data2[which(data2$stability=="stable" & data2$copy=="low"),x]))
	m5<-median(as.matrix(data2[-(which(data2$stability=="stable" & data2$copy=="low")),x]))
	m6<-wilcox.test(max(m1,m2), min(m1,m2), paired = TRUE, alternative = "greater")$p.value

	val<-c(colnames(data2)[x],m1,m2, m3,m4,m5,m6)
	data.mat <- rbind(data.mat, val)
}
prmatrix(data.mat,quote=FALSE, rowlab=rep("",8), collab=rep("",7))

data.mat <- c("Field","stable.mean","unstable.mean","p.val.mean","stable.median","unstable.median","p.val.median")
for (x in 7:length(colnames(data))){
	m1<-mean(as.matrix(data2[which(data2$stability=="stable"),x]))
	m2<-mean(as.matrix(data2[-(which(data2$stability=="stable")),x]))
	m3<-wilcox.test(max(m1,m2), min(m1,m2), paired = TRUE, alternative = "greater")$p.value

	m4<-median(as.matrix(data2[which(data2$stability=="stable"),x]))
	m5<-median(as.matrix(data2[-(which(data2$stability=="stable")),x]))
	m6<-wilcox.test(max(m1,m2), min(m1,m2), paired = TRUE, alternative = "greater")$p.value

	val<-c(colnames(data2)[x],m1,m2, m3,m4,m5,m6)
	data.mat <- rbind(data.mat, val)
}
prmatrix(data.mat,quote=FALSE, rowlab=rep("",8), collab=rep("",7))


data2["phenotype"]<-rep('1', times=44)
data2$phenotype[data2$Favorability.numeric!=3]<-0

data2["stab"]<-rep('1', times=44)
data2$stab[data2$Favorability.numeric<0]<-0

data3<-read.table("janssen_complete_analysis.txt", header=TRUE, sep="\t")
data3["phenotype"]<-rep('1', times=24)
data3$phenotype[data3$Favorability.numeric!=3]<-0
data3$phenotype<-as.numeric(data3$phenotype)

data3["stab"]<-rep('1', times=24)
data3$stab[data3$Favorability.numeric<0]<-0
data3$stab<-as.numeric(data3$stab)
data3<-transform(data3, min.expression.dis= pmin(High.expression.bin.to.the.left, High.expression.bin.to.the.right))
data3<-transform(data3, min.variant.dis= pmin(High.variant.bin.to.the.left, High.variant.bin.to.the.right))

mylogit <- glm(stab ~ Variant.Profile + min.expression.dis + state +Structural.Variation +  TG.TG.Fusions, data = data3, family = "binomial" )

x<-data3[which(data3$stability=="stable"),26]
y<-data3[which(data3$stability=="unstable"),26]
mean(x)
mean(y)
t.test(y, x,  alternative = "greater")





data<-read.table("complete_analysis.txt", header=TRUE, sep="\t")
data["phenotype"]<-rep('1', times=44)
data$phenotype[data$Favorability.numeric!=3]<-0
data$phenotype<-as.numeric(data$phenotype)

data["stab"]<-rep('1', times=44)
data$stab[data$Favorability.numeric<0]<-0
data$stab<-as.numeric(data$stab)
data<-transform(data, min.expression.dis= pmin(High.expression.bin.to.the.left, High.expression.bin.to.the.right))
data<-transform(data, min.variant.dis= pmin(High.variant.bin.to.the.left, High.variant.bin.to.the.right))
data<-transform(data, min.reg.dis= pmin(Regulatory.region.left, Regulatory.region.right))

mylogit <- glm(stab ~ Variant.Profile + min.expression.dis + state + Structural.Variation +  TG.TG.Fusions, data = data, family = "binomial" )
mylogit <- glm(stab ~ Variant.Profile + min.expression.dis + state, data = data, family = "binomial")


x<-data[which(data$stability=="stable"),22]
y<-data[which(data$stability=="unstable"),22]
mean(x)
mean(y)
wilcox.test(y, x, paired = TRUE, alternative = "greater")


data.bad<-read.table("bad_features.txt", header=TRUE)
data.bad.m<-melt(data.bad, id.vars=c("cell.line", "copy.number", "stability", "favorability", "favorability.numeric"))
library(reshape2)
p<-ggplot(data.bad.m, aes(x=cell.line, y=value, fill=factor(variable))) + 
geom_bar(stat="identity", position="dodge") + 
theme_light() + 
theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
labs(x = "Cell lines", fill = "Factors") + 
scale_x_discrete(limits=data.bad$cell.line) + 
ggsave("bad_features.pdf", p)

summary(data[which(data$stab==1),"min.reg.dis"])
summary(data[which(data$stab==0),"min.reg.dis"])
summary(data[which(data$stab==1),"min.expression.dis"])
summary(data[which(data$stab==0),"min.expression.dis"])
summary(data[which(data$stab==1),"min.variant.dis"])
summary(data[which(data$stab==0),"min.variant.dis"])

boxplot(data[which(data$stab==1),"min.expression.dis"], data[which(data$stab==0),"min.expression.dis"])
boxplot(data[which(data$stab==1),"min.variant.dis"], data[which(data$stab==0),"min.variant.dis"])
boxplot(data[which(data$stab==1),"min.reg.dis"], data[which(data$stab==0),"min.reg.dis"])

boxplot(data[which(data$phenotype==1),"min.expression.dis"], data[which(data$phenotype==0),"min.expression.dis"])
boxplot(data[which(data$phenotype==1),"min.variant.dis"], data[which(data$phenotype==0),"min.variant.dis"])
boxplot(data[which(data$phenotype==1),"min.reg.dis"], data[which(data$phenotype==0),"min.reg.dis"])

summary(data[which(data$phenotype==1),"min.reg.dis"])
summary(data[which(data$phenotype==0),"min.reg.dis"])
summary(data[which(data$phenotype==1),"min.expression.dis"])
summary(data[which(data$phenotype==0),"min.expression.dis"])
summary(data[which(data$phenotype==1),"min.variant.dis"])
summary(data[which(data$phenotype==0),"min.variant.dis"])

p<-ggplot(data, aes(x=factor(stab), y=min.expression.dis)) + geom_violin(aes(fill=factor(stab)), trim = FALSE)
p<-p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)  + stat_summary(fun.y=median, geom="point", size=2, color="red") + labs(x = "Stability", fill = "Stability") 
ggsave("exp_dis_vioplot.pdf", p)
p2<-ggplot(data, aes(x=factor(stab), y=min.expression.dis)) + geom_boxplot(aes(fill=factor(stab))) + labs(x = "Stability", fill = "Stability") + theme_light()
ggsave("exp_dis_boxplot.pdf", p2)


a<-ggplot(data, aes(x=factor(stab), y=Variant.Profile)) + geom_boxplot()+ xlab("Stability") + ylab("Variant profile at the site") + theme_light()
b<-ggplot(data, aes(x=factor(stab), y=Expression.Profile)) + geom_boxplot()+ xlab("Stability") + ylab("Expression profile at the site")+ theme_light()
c<-ggplot(data, aes(x=factor(stab), y=state)) + geom_boxplot()+ xlab("Stability") + ylab("Chromatin state at/around the site")+ theme_light()
d<-ggplot(data, aes(x=factor(stab), y=min.variant.dis)) + geom_boxplot()+ xlab("Stability") + ylab("Min. distance from high variant bin") + theme_light()
e<-ggplot(data, aes(x=factor(stab), y=min.expression.dis)) + geom_boxplot()+ xlab("Stability") + ylab("Min. distance from high expression bin") + theme_light()
f<-ggplot(data, aes(x=factor(stab), y=min.reg.dis)) + geom_boxplot()+ xlab("Stability") + ylab("Min. distance from regulatory region") + theme_light()
p3<-grid.arrange(a,b,c,d,e,f, nrow = 2)
ggsave("factors_boxplot.pdf", p3)


a<-ggplot(data, aes(x=factor(phenotype), y=Variant.Profile)) + geom_boxplot()+ xlab("Phenotype") + ylab("Variant profile at the site") + theme_light()
b<-ggplot(data, aes(x=factor(phenotype), y=Expression.Profile)) + geom_boxplot()+ xlab("Phenotype") + ylab("Expression profile at the site")+ theme_light()
c<-ggplot(data, aes(x=factor(phenotype), y=state)) + geom_boxplot()+ xlab("Phenotype") + ylab("Chromatin state at/around the site")+ theme_light()
d<-ggplot(data, aes(x=factor(phenotype), y=min.variant.dis)) + geom_boxplot()+ xlab("Phenotype") + ylab("Min. distance from high variant bin") + theme_light()
e<-ggplot(data, aes(x=factor(phenotype), y=min.expression.dis)) + geom_boxplot()+ xlab("Phenotype") + ylab("Min. distance from high expression bin") + theme_light()
f<-ggplot(data, aes(x=factor(phenotype), y=min.reg.dis)) + geom_boxplot()+ xlab("Phenotype") + ylab("Min. distance from regulatory region") + theme_light()
p3<-grid.arrange(a,b,c,d,e,f, nrow = 2)
ggsave("factors_phenotype_boxplot.pdf", p3)


# p1<-ggplot(data, aes(x=high.variant.bin.upstream, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + ggtitle("High variant bin to the left")  + theme(plot.title = element_text(size=8)) + theme(text = element_text(size=6)) + theme_light() + scale_shape_manual(values=seq(0,15))
# p2<-ggplot(data, aes(x=high.variant.bin.downstream, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + ggtitle("High variant bin to the right")  + theme(plot.title = element_text(size=8)) + theme(text = element_text(size=6)) + theme_light() + scale_shape_manual(values=seq(0,15))
# p3<-ggplot(data, aes(x=high.expression.bin.upstream, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + ggtitle("High expression bin to the left")  + theme(plot.title = element_text(size=8)) + theme(text = element_text(size=6)) + theme_light() + scale_shape_manual(values=seq(0,15))
# p4<-ggplot(data, aes(x=high.expression.bin.downstream, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + ggtitle("High expression bin to the right")  + theme(plot.title = element_text(size=8)) + theme(text = element_text(size=6)) + theme_light() + scale_shape_manual(values=seq(0,15))
# p5<-ggplot(data, aes(x=regulatory.region.left, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + ggtitle("High expression bin to the left")  + theme(plot.title = element_text(size=8)) + theme(text = element_text(size=6)) + theme_light() + scale_shape_manual(values=seq(0,15))
# p6<-ggplot(data, aes(x=regulatory.region.right, y=Favorability.numeric,  color=Favorability)) + geom_point(size=2, position=position_jitter(h=0.15,w=0.15), show.legend=F) + ggtitle("High expression bin to the right")  + theme(plot.title = element_text(size=8)) + theme(text = element_text(size=6)) + theme_light() + scale_shape_manual(values=seq(0,15))
# grid.arrange(p1, p2, p3, p4,p5,p6, nrow = 3)


# m1<-mean(as.matrix(data[which(data$stability=="stable" & data$copy=="low"),11:12]))
# m2<-mean(as.matrix(data[-(which(data$stability=="stable" & data$copy=="low")),11:12]))
# wilcox.test(m1, m2, paired = TRUE, alternative = "greater")
# m1
# m2
# m1<-mean(as.matrix(data[which(data$stability=="stable" & data$copy=="low"),13:14]))
# m2<-mean(as.matrix(data[-(which(data$stability=="stable" & data$copy=="low")),13:14]))
# wilcox.test(m1, m2, paired = TRUE, alternative = "greater")
# m1
# m2
# m1<-mean(as.matrix(data[which(data$stability=="stable" & data$copy=="low"),c(16,18)]))
# m2<-mean(as.matrix(data[-(which(data$stability=="stable" & data$copy=="low")),c(16,18)]))
# wilcox.test(m2, m1, paired = TRUE, alternative = "greater")
# m1
# m2
# m1<-median(as.matrix(data[which(data$stability=="stable" & data$copy=="low"),c(16,18)]))
# m2<-median(as.matrix(data[-(which(data$stability=="stable" & data$copy=="low")),c(16,18)]))
# wilcox.test(m1, m2, paired = TRUE, alternative = "greater")
# m1
# m2


# m3=mean(as.matrix(data[which(data$stability=="stable" & data$copy=="low"),11:12]))
# m4=mean(as.matrix(data[-(which(data$stability=="stable" & data$copy=="low")),11:12]))
# m5=mean(as.matrix(data[which(data$stability=="stable" & data$copy=="low"),13:14]))
# m6=mean(as.matrix(data[-(which(data$stability=="stable" & data$copy=="low")),13:14]))
# m7=mean(as.matrix(data[which(data$stability=="stable" & data$copy=="low"),c(16,18)]))
# m8=mean(as.matrix(data[-(which(data$stability=="stable" & data$copy=="low")),c(16,18)]))
# m<-c(m3,m4,m5,m6,m7,m8)
# m
# 57.70770 37.14706 
# 348.21913 483.88212 
# 27.67556 30.14318
# m3=mean(as.matrix(data[which(data$stability=="stable"),11:12]))
# m4=mean(as.matrix(data[-(which(data$stability=="stable")),11:12]))
# m5=mean(as.matrix(data[which(data$stability=="stable"),13:14]))
# m6=mean(as.matrix(data[-(which(data$stability=="stable")),13:14]))
# m7=mean(as.matrix(data[which(data$stability=="stable"),c(16,18)]))
# m8=mean(as.matrix(data[-(which(data$stability=="stable")),c(16,18)]))
# m<-c(m3,m4,m5,m6,m7,m8)
# m
# 50.51842 47.50000 
# 334.12847 600.15209 
# 27.26961 32.70700
