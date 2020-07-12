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
