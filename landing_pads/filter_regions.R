require("GenomicRanges")
fav<-read.table("results2/filter_regions_v1_all_parsed_filtered2.txt", header=F, sep="\t")
colnames(fav)<-c("chr", "start", "stop", "high.exp", "high.var", "qui", "rep", "enh", "prom", "act", "anno", "chr.ass", "size")

site_fav<-read.table("src2/favorable_ncbi.txt", header=F, sep="\t")
site_unfav<-read.table("src2/unfavorable_ncbi.txt", header=F, sep="\t")

sites_old_sel<-read.table("sel_fav_ncbi.txt", header=T)
sites_old_sel_unfav<-read.table("sel_unfav_ncbi.txt", header=T)
colnames(site_fav)<-c("chr", "start", "stop", "cell line", "stability")
colnames(site_unfav)<-c("chr", "start", "stop", "cell line", "stability")

fav$chr.ass<-factor(fav$chr.ass, levels=c("2","X","7","3","5","1","4","8","6","10","9","unplaced"))

### Assign max value to the parameters on scaffolds without any negative parameter
fav[is.na(fav$qui),"qui"]<-max(fav[!is.na(fav$qui),"qui"])
fav[is.na(fav$rep),"rep"]<-max(fav[!is.na(fav$rep),"rep"])
fav[is.na(fav$high.var),"high.var"]<-max(fav[!is.na(fav$high.var),"high.var"])

### Measure quantiles or scale by max value
fav.q<-fav
fav.q$high.exp <- fav.q$high.exp/max(fav.q$high.exp)
fav.q$high.var <- fav.q$high.var/max(fav.q$high.var)
fav.q$qui <- fav.q$qui/max(fav.q$qui)
fav.q$rep <- fav.q$rep/max(fav.q$rep)
fav.q$enh <- fav.q$enh/max(fav.q$enh)
fav.q$prom <- fav.q$prom/max(fav.q$prom)
fav.q$enh <- fav.q$enh/max(fav.q$enh)
fav.q$prom <- fav.q$prom/max(fav.q$prom)
fav.q$act <- fav.q$act/max(fav.q$act)
fav.q$size <- fav.q$size/max(fav.q$size)

### Create score 
fav.q$score <- fav.q$size*(fav.q$high.var + fav.q$rep + fav.q$size + (4 - (fav.q$high.exp + fav.q$enh + fav.q$prom + fav.q$act)))

### Sort the filtered regions
fav.q.sort <- fav.q[with(fav.q, order(-score,chr.ass,anno)), ]
fav.q.sort$rank <- c(1:dim(fav.q.sort)[1])
fav.q.sort$loci<-paste(paste(fav.q.sort$chr, fav.q.sort$start,sep=":"),fav.q.sort$stop, sep="..")
#### Write to table
write.table(fav.q.sort, "results2/favorable_sorted.txt",sep="\t",quote=F,row.names=F,col.names=T)

fav$loci<-paste(paste(fav$chr, fav$start,sep=":"),fav$stop, sep="..")
rownames(fav)<-fav$loci
rownames(fav.q.sort)<-fav.q.sort$loci
final.fav<-merge(fav, fav.q.sort, by="row.names")
final.fav<-check[with(check, order(rank)),]
write.table(final.fav, "results2/favorable_sorted_final.txt",sep="\t",quote=F,row.names=F,col.names=T)



head(fav.q.sort, n=50L)
fav.gr <- makeGRangesFromDataFrame(fav, keep.extra.columns=TRUE)
fav.q.sort.gr <- makeGRangesFromDataFrame(fav.q.sort, keep.extra.columns=TRUE)
site_fav.gr <- makeGRangesFromDataFrame(site_fav, keep.extra.columns=TRUE)
site_unfav.gr <- makeGRangesFromDataFrame(site_unfav, keep.extra.columns=TRUE)
sites_old_sel.gr <- makeGRangesFromDataFrame(sites_old_sel, keep.extra.columns=TRUE)

### Overlap with known favorable sites
site_fav.ol<-findOverlaps(site_fav.gr, fav.q.sort.gr)
site_fav.ol_info<-fav.q.sort.gr[subjectHits(site_fav.ol),]
site_fav.ol_info
site_fav.gr[queryHits(site_fav.ol),]

ori.ol<-findOverlaps(site_fav.ol_info, fav.gr)
fav.gr[subjectHits(ori.ol),]

### Overlap with known unfavorable sites
site_unfav.ol<-findOverlaps(site_unfav.gr, fav.q.sort.gr)
fav.q.sort.gr[subjectHits(site_unfav.ol),]

### Overlap with old selected favorable sites
site_old.ol<-findOverlaps(sites_old_sel.gr, fav.q.sort.gr)
fav.q.sort.gr[subjectHits(site_old.ol),]

############################################################

unfav<-read.table("results2/filter_regions_v2_all_parsed_filtered2.txt", header=F, sep="\t")
colnames(unfav)<-c("chr", "start", "stop", "high.exp", "high.var", "qui", "rep", "enh", "prom", "act", "anno", "chr.ass", "size")
unfav$chr.ass<-factor(unfav$chr.ass, levels=c("2","X","7","3","5","1","4","8","6","10","9","unplaced"))

### Assign max value to the parameters on scaffolds without any negative parameter
unfav[is.na(unfav$qui),"qui"]<-max(unfav[!is.na(unfav$qui),"qui"])
unfav[is.na(unfav$rep),"rep"]<-max(unfav[!is.na(unfav$rep),"rep"])
unfav[is.na(unfav$high.var),"high.var"]<-max(unfav[!is.na(unfav$high.var),"high.var"])

### Measure quantiles or scale by max value
unfav.q<-unfav
unfav.q$high.exp <- unfav.q$high.exp/max(unfav.q$high.exp)
unfav.q$high.var <- unfav.q$high.var/max(unfav.q$high.var)
unfav.q$qui <- unfav.q$qui/max(unfav.q$qui)
unfav.q$rep <- unfav.q$rep/max(unfav.q$rep)
unfav.q$enh <- unfav.q$enh/max(unfav.q$enh)
unfav.q$prom <- unfav.q$prom/max(unfav.q$prom)
unfav.q$enh <- unfav.q$enh/max(unfav.q$enh)
unfav.q$prom <- unfav.q$prom/max(unfav.q$prom)
unfav.q$act <- unfav.q$act/max(unfav.q$act)
unfav.q$size <- unfav.q$size/max(unfav.q$size)

### Create score 
unfav.q$score <- unfav.q$high.var + unfav.q$rep + unfav.q$size + (4 - (unfav.q$high.exp + unfav.q$enh + unfav.q$prom + unfav.q$act))

### Sort the filtered regions
unfav.q.sort <- unfav.q[with(unfav.q, order(score,rev(chr.ass),rev(anno))), ]
unfav.q.sort$rank <- c(1:dim(unfav.q.sort)[1])
unfav.q.sort$loci<-paste(paste(unfav.q.sort$chr, unfav.q.sort$start,sep=":"), unfav.q.sort$stop, sep="..")

#### Write to table
write.table(unfav.q.sort, "results2/unfavorable_sorted.txt",sep="\t",quote=F,row.names=F,col.names=T)

unfav$loci<-paste(paste(unfav$chr, unfav$start,sep=":"),unfav$stop, sep="..")
rownames(unfav)<-unfav$loci
rownames(unfav.q.sort)<-unfav.q.sort$loci
final.unfav<-merge(unfav, unfav.q.sort, by="row.names")
final.unfav<-check[with(check, order(rank)),]
write.table(final.unfav, "results2/unfavorable_sorted_final.txt",sep="\t",quote=F,row.names=F,col.names=T)



head(unfav.q.sort, n=50L)

unfav.gr <- makeGRangesFromDataFrame(unfav, keep.extra.columns=TRUE)
unfav.q.sort.gr <- makeGRangesFromDataFrame(unfav.q.sort, keep.extra.columns=TRUE)
site_fav.gr <- makeGRangesFromDataFrame(site_fav, keep.extra.columns=TRUE)
site_unfav.gr <- makeGRangesFromDataFrame(site_unfav, keep.extra.columns=TRUE)
sites_old_sel_unfav.gr <- makeGRangesFromDataFrame(sites_old_sel_unfav, keep.extra.columns=TRUE)

### Overlap with known unfavorable sites
site_unfav.ol<-findOverlaps(site_unfav.gr, unfav.q.sort.gr)
unfav.q.sort.gr[subjectHits(site_unfav.ol),]
site_unfav.gr[queryHits(site_unfav.ol),]

### Overlap with known favorable sites
site_unfav_fav.ol<-findOverlaps(site_fav.gr, unfav.q.sort.gr)
unfav.q.sort.gr[subjectHits(site_unfav_fav.ol),]
site_fav.gr[queryHits(site_unfav_fav.ol),]

### Overlap with old selected unfavorable sites
site_old_unfav.ol<-findOverlaps(sites_old_sel_unfav.gr, unfav.q.sort.gr)
unfav.q.sort.gr[subjectHits(site_old_unfav.ol),]
sites_old_sel_unfav.gr[queryHits(site_old_unfav.ol),]

