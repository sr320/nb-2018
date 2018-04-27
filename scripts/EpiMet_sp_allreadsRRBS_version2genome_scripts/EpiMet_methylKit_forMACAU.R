library(methylKit)
source("http://www.bioconductor.org/biocLite.R")
library(data.table)
library(GenomicRanges)
library(mixtools)
library(data.table)
library(rtracklayer)
library(pheatmap)
setwd("~/clusteroutput/EpiMet_Sp_version2genome_released/")
file.list <- list("EpiMet1_sp_allreads_reformcov.txt",
                  "EpiMet2_sp_allreads_reformcov.txt",
                  "EpiMet3_sp_allreads_reformcov.txt",
                  "EpiMet4_sp_allreads_reformcov.txt",
                  "EpiMet5_sp_allreads_reformcov.txt",
                  "EpiMet6_sp_allreads_reformcov.txt",
                  "EpiMet9_sp_allreads_reformcov.txt",
                  "EpiMet10_sp_allreads_reformcov.txt",
                  "EpiMet11_sp_allreads_reformcov.txt",
                  "EpiMet12_sp_allreads_reformcov.txt",
                  "EpiMet13_sp_allreads_reformcov.txt",
                  "EpiMet14_sp_allreads_reformcov.txt",
                  "EpiMet15_sp_allreads_reformcov.txt",
                  "EpiMet17_sp_allreads_reformcov.txt",
                  "EpiMet18_sp_allreads_reformcov.txt",
                  "EpiMet19_sp_allreads_reformcov.txt",
                  "EpiMet20_sp_allreads_reformcov.txt",
                  "EpiMet21_sp_allreads_reformcov.txt",
                  "EpiMet24_sp_allreads_reformcov.txt",
                  "EpiMet25_sp_allreads_reformcov.txt",
                  "EpiMet26_sp_allreads_reformcov.txt",
                  "EpiMet27_sp_allreads_reformcov.txt",
                  "EpiMet30_sp_allreads_reformcov.txt",
                  "EpiMet31_sp_allreads_reformcov.txt",
                  "EpiMet32_sp_allreads_reformcov.txt",
                  "EpiMet33_sp_allreads_reformcov.txt",
                  "EpiMet34_sp_allreads_reformcov.txt",
                  "EpiMet35_sp_allreads_reformcov.txt",
                  "EpiMet36_sp_allreads_reformcov.txt",
                  "EpiMet37_sp_allreads_reformcov.txt",
                  "EpiMet38_sp_allreads_reformcov.txt",
                  "EpiMet39_sp_allreads_reformcov.txt",
                  "EpiMet40_sp_allreads_reformcov.txt",
                  "EpiMet41_sp_allreads_reformcov.txt",
                  "EpiMet42_sp_allreads_reformcov.txt",
                  "EpiMet43_sp_allreads_reformcov.txt",
                  "EpiMet44_sp_allreads_reformcov.txt",
                  "EpiMet45_sp_allreads_reformcov.txt",
                  "EpiMet46_sp_allreads_reformcov.txt",
                  "EpiMet47_sp_allreads_reformcov.txt",
                  "EpiMet48_sp_allreads_reformcov.txt",
                  "EpiMet49_sp_allreads_reformcov.txt",
                  "EpiMet50_sp_allreads_reformcov.txt",
                  "EpiMet51_sp_allreads_reformcov.txt",
                  "EpiMet52_sp_allreads_reformcov.txt",
                  "EpiMet53_sp_allreads_reformcov.txt",
                  "EpiMet54_sp_allreads_reformcov.txt",
                  "EpiMet55_sp_allreads_reformcov.txt",
                  "EpiMet56_sp_allreads_reformcov.txt",
                  "EpiMet57_sp_allreads_reformcov.txt",
                  "EpiMet58_sp_allreads_reformcov.txt",
                  "EpiMet59_sp_allreads_reformcov.txt",
                  "EpiMet60_sp_allreads_reformcov.txt",
                  "EpiMet61_sp_allreads_reformcov.txt",
                  "EpiMet62_sp_allreads_reformcov.txt",
                  "EpiMet63_sp_allreads_reformcov.txt",
                  "EpiMet64_sp_allreads_reformcov.txt",
                  "EpiMet65_sp_allreads_reformcov.txt",
                  "EpiMet66_sp_allreads_reformcov.txt",
                  "EpiMet67_sp_allreads_reformcov.txt") 
myobj<-methRead( file.list,pipeline=list(fraction=TRUE,chr.col=2,start.col=3,end.col=3,
                                     coverage.col=5,strand.col=4,freqC.col=6 ),
             sample.id=list("T_n_64",
                            "T_k_73",
                            "T_l_81",
                            "S_d_93",
                            "S_l_95",
                            "S_i_115",
                            "T_d_83",
                            "T_i_75",
                            "T_d_84",
                            "T_f_85",
                            "S_c_107",
                            "S_k_109",
                            "S_d_111",
                            "T_d_61",
                            "T_c_89",
                            "S_n_104",
                            "S_f_113",
                            "S_d_121",
                            "T_n_82",
                            "T_p_65",
                            "S_n_97",
                            "S_p_122",
                            "T_i_63",
                            "S_a_119",
                            "T_h_78",
                            "T_m_80",
                            "T_h_76",
                            "T_r_79",
                            "S_n'_123",
                            "T_i_71",
                            "T_f_90",
                            "S_a_92",
                            "S_c_101",
                            "S_c_102",
                            "T_f_68",
                            "S_o'_116",
                            "T_i_86",
                            "T_r_72",
                            "T_m_87",
                            "T_b_77",
                            "S_g_98",
                            "S_o'_114",
                            "T_i_74",
                            "S_l_112",
                            "S_c_99",
                            "S_q_118",
                            "S_k_91",
                            "S_d_117",
                            "S_c_106",
                            "S_n_108",
                            "T_h_66",
                            "S_q_96",
                            "S_a'_94",
                            "S_c_100",
                            "S_l_110",
                            "T_f_62",
                            "T_m_67",
                            "T_m_88",
                            "T_s_70",
                            "T_o_69"),assembly="version2_released",
             treatment=c(1,1,1,0,0,0,1,1,1,1,0,0,0,1,1,0,0,0,1,1,0,0,1,0,1,1,1,1,0,1,1,0,0,0,1,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,1,1))


#NOTE: amanda says she uses those with 50% of data, methylkit will let me only do it by group, so will say 50% per group n=30 per group, so 15)
meth<-unite(myobj, min.per.group=15L)
nrow(meth)
head(meth)


#generate a df from the methylobj
meth_df<-data.frame(meth)
head(meth_df)


# Sample ids from the meth data obj (used later)
sample.ids <- slot(meth, "sample.ids")

#add a column with positionID
meth_df$site<-paste(meth_df$chr,meth_df$start,meth_df$start,sep=".")
head(meth_df)
#want to add a filter here ... only retain sites where median read depth is >10

tmp.coverage <- meth_df[, grep(pattern="coverage", x=colnames(meth_df))]
tmp.coverage.medians <- apply(X=tmp.coverage, FUN=median, MARGIN=1, na.rm=T)
#tmp.coverage.means <- apply(X=tmp.coverage, FUN=mean, MARGIN=1, na.rm=T)

# Meth_df filtered by medians across all samples
meth_df.filtered_medians <- meth_df[tmp.coverage.medians > 10, ]

#get percent methylation per site
percents<-percMethylation(meth, rowids=T)
head(percents)
write.table(percents,"percents_sp_allreads_60indv_15indvmin.txt", sep = "\t", row.names = T, col.names = T, quote=F)
#want to make a white list of sites where 1) filter those sites that fall in bottom 5% of dataset in terms of variance, 2) mean methylation >90%, 3) mean methylation <10%

# Variance Whitelist
tmp.percent.variances <- apply(X=percents, FUN=var, MARGIN=1, na.rm=T)
percent.variance.cutoff <- quantile(tmp.percent.variances, probs=0.05)
whitelist.variances <- rownames(percents)[tmp.percent.variances > percent.variance.cutoff]

# Hyper Methylated Whitelist
tmp.percent.means <- apply(X=percents, FUN=mean, MARGIN=1, na.rm=T)
percent.hyper.cutoff <- 90
whitelist.hyper <- rownames(percents)[tmp.percent.means < percent.hyper.cutoff]

# Hypo Methylated Whitelist
percent.hypo.cutoff <- 10
whitelist.hypo <- rownames(percents)[tmp.percent.means > percent.hypo.cutoff]

# Union whitelists
whitelist.percents.final <- intersect(intersect(whitelist.variances, whitelist.hyper), whitelist.hypo)

#####working out how many sites are lost PER filter########
whitelist.after.low.variance.removed<-#1419554
whitelist.after.hypo.removed<-intersect(whitelist.variances,whitelist.hypo) #1359452
whitelist.after.hyper.removed<-intersect(whitelist.after.hypo.removed,whitelist.hyper) #84661

#join the whitelist to the coverage values to only retain those sites which pass the filters
meth_df.filtered.final <- meth_df.filtered_medians[meth_df.filtered_medians$site %in% whitelist.percents.final, ]

# make 2 new tables from the coverage values - 1 with coverage and 1 with methylation values
coverage.table <- meth_df.filtered.final[, c(grep(pattern="site", x=colnames(meth_df.filtered.final)), 
                                             grep(pattern="coverage", x=colnames(meth_df.filtered.final)))]
colnames(coverage.table) <- c("site", sample.ids)

# numC Table
numC.table <- meth_df.filtered.final[, c(grep(pattern="site", x=colnames(meth_df.filtered.final)), 
                                         grep(pattern="numCs", x=colnames(meth_df.filtered.final)))]
colnames(numC.table) <- c("site", sample.ids)

#write the 2 files as .txt
write.table(coverage.table,"totalCounts_sp_allreads60indv.txt", sep = "\t", row.names = F, col.names = T, quote=F)
write.table(numC.table,"methCounts_sp_allreads60indv.txt", sep = "\t", row.names = F, col.names = T, quote=F)




####after MACAU is performed, read in the output table to do multiple testing correction
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)
setwd("~/clusteroutput/EpiMet_Sp_version2genome_released/MACAUfiles/output_EpiMet_sp_allreads/")
macau_out<-read.table("output_EpiMet_sp_allreads60indv.assoc.txt", sep = '\t',header=T)
head(macau_out)
nrow(macau_out)
#how many CpG sites have pvalue < 0.01?
sum(macau_out$pvalue <= 0.01499999)
macau_out$fdr <- lfdr(macau_out$pvalue)
#how many CpG sites have FDR < 0.1?
sum(macau_out$fdr <= 0.1499999)
write.table(macau_out,"output_sp_allreads60indv", sep = "\t", row.names = F, col.names = T, quote=F)

#hist(macau_out$pvalue)
#hist(macau_out$h)
##############################################Downstream Stuff######################################################
###making a heatmap of the significant CpGs
setwd("~/clusteroutput/EpiMet_Sp_version2genome_released/MACAUfiles/output_EpiMet_sp_allreads/")
sig_sites<-read.table("sig_sitesFDR.05_sp_allreads60indv.txt", sep = '\t',header=T)
head(sig_sites)
#You can then do something like this
row.names(percents) %in% sig_sites
#This will give you a list of T/F values which you can then subset your table.  For example
subset_percents<-percents[row.names(percents) %in% sig_sites$sig_sites, ]
subset_percents
write.table(subset_percents,"MACAUoutput_notrim60indv_FDR_sigsites_percents", sep = "\t", row.names = F, col.names = T, quote=F)
####making a heatmap of methylkit sig sites
library(gplots)
library(RColorBrewer)
subset_MACAU<-heatmap.2(as.matrix(subset_percents),na.color="grey",col=brewer.pal(9,"Reds"),trace=c("none"))

###this is methylkit analysis
#You can then do something like this
row.names(percents) %in% myDiffq_df
#This will give you a list of T/F values which you can then subset your table.  For example

subset_percents<-percents[row.names(percents) %in% myDiffq_df$site, ]
head(subset_percents)
library(gplots)
library(RColorBrewer)
subset_methylkit<-heatmap.2(as.matrix(subset_percents),na.color="grey",col=brewer.pal(9,"Reds"),trace=c("none"))

#significant sites from MACAU additional breakdown
MACAUsigPerc<-read.table("MACAUoutput_notrim60indv_FDR_sigsites_percents", sep = "\t", header = T)
head(MACAUsigPerc)
one<-ggplot(MACAUsigPerc, aes(x=family, y=site1)) + geom_boxplot()
two<-ggplot(MACAUsigPerc, aes(x=family, y=site2)) + geom_boxplot()
three<-ggplot(MACAUsigPerc, aes(x=family, y=site3)) + geom_boxplot()
four<-ggplot(MACAUsigPerc, aes(x=family, y=site4)) + geom_boxplot()
five<-ggplot(MACAUsigPerc, aes(x=family, y=site5)) + geom_boxplot()
six<-ggplot(MACAUsigPerc, aes(x=family, y=site6)) + geom_boxplot()
seven<-ggplot(MACAUsigPerc, aes(x=family, y=site7)) + geom_boxplot()
require(cowplot)
plot_grid(one,two,three,four,five,six,seven, align = 'v', nrow = 7, ncol = 1,hjust=-1.5,label_size=20)



########methylkit analysis########
#########ANALYSIS#########################
###########TILING#########################
# Function
tileit <- function(object, win.size=100, step.size=100, cov.bases=20) {
  g.meth <- as(object,"GRanges")
  chrs <- as.character(unique(seqnames(g.meth)))
  grl <- split(g.meth, seqnames(g.meth))
  max.length <- max(end(grl))
  numTiles <- floor((max.length - (win.size - step.size))/step.size) + 1
  starts <- unlist(sapply(numTiles,function(x) 1+0:(x-1)*step.size))
  ranges <- IRanges(start=starts, width=rep(win.size,length(starts)))
  all.wins <- GRanges(seqnames=Rle(chrs, numTiles), ranges=ranges)
  rcounts <- regionCounts(object, all.wins, 0, strand.aware=FALSE)
  rcounts.filtered <- rcounts[rcounts$coverage >= cov.bases, ]
  return(rcounts.filtered)
}

# Code to run function on methylKitList
new.list <- lapply(myobj, tileit)
myobj.tiled <- new("methylRawList", new.list, treatment=myobj@treatment)
head(myobj.tiled)



#normalize and unite
meth<-unite(myobj.tiled,min.per.group=30L)
nrow(meth)
head(meth)

#getpercents
percents<-percMethylation(meth, rowids=T)
head(percents)
write.table(percents,"percents_100bpTiles_sp_allreads60indv.txt", sep = "\t", row.names = T, col.names = F, quote=F)


###Correlation####
#check correlation between samples
#careful if you run this with a lot of samples because it takes FOREVER
#getCorrelation(meth,plot=T)
#this makes same graph, but saves it higher resolution!
#pdf(file="individual_graph.pdf")
#getCorrelation(meth,plot=T)
#dev.off()

####clustering#####
#dendrogram
#lots of options for distance and method changes, 'correlation' and 'ward' are defaults
clusterSamples(meth, dist="correlation", method="average", plot=TRUE)
#nrow(meth)
#?clusterSamples
#PCA
#need to check these settings..
?PCASamples
PCASamples(meth, screeplot=TRUE)
PCASamples(meth,adj.lim=c(0.4,0.2))
PCASamples(meth,screeplot=FALSE, adj.lim=c(0.4,0.2),
           scale=TRUE,center=TRUE,comp=c(3,4),transpose=TRUE,sd.filter=TRUE,
           sd.threshold=0.7,filterByQuantile=TRUE,obj.return=FALSE)

#testing for batch effects
batch<-read.table("covariates_sequencer_totalmeth_family.txt",header=F,sep = '\t')
head(batch)
as<-assocComp(mBase=meth,batch)
as
newObj=removeComp(meth,comp=1)
newObj2=removeComp(newObj,comp = 1)
newObj3=removeComp(newObj2,comp=1)
PCASamples(newObj2,screeplot=FALSE, adj.lim=c(0.4,0.2),
           scale=TRUE,center=TRUE,comp=c(1,2),transpose=TRUE,sd.filter=TRUE,
           sd.threshold=0.7,filterByQuantile=TRUE,obj.return=FALSE)

#calculate diffs
?calculateDiffMeth
covariate<-read.table("dam_sire_covar.txt", header=T,sep = '\t')
head(covariate)
covariate <- as.data.frame(covariate)
head(covariate)


myDiff<-calculateDiffMeth(meth,num.cores = 8,slim=TRUE, covariates = covariate)

head(myDiff)
hist(myDiff$diff.meth)
#sig diffs
myDiffq <- getMethylDiff(myDiff, difference = 25, qvalue = 0.01)
write.table(myDiffq,"methylKit_tiles_EpiMet_DamSireCovar_20percq0.01.txt", sep = "\t", row.names = T, col.names = F, quote=F)

nrow(myDiffq)
head(myDiffq)

sig_sites<-read.table("MACAUfiles/sigsites_EpiMet_Sp_MACAUpval0.001.txt", sep = '\t',header=F)
head(sig_sites)
#You can then do something like this
row.names(percents)
row.names(percents) %in% sig_sites
#This will give you a list of T/F values which you can then subset your table.  For example
subset_percents<-percents[row.names(percents) %in% sig_sites$V1, ]
subset_percents
write.table(subset_percents,"methylKit_tiles_EpiMet_DamSireCovar_20percq0.01.percents.txt", sep = "\t", row.names = F, col.names = T, quote=F)
####making a heatmap of methylkit sig sites
library(gplots)
library(RColorBrewer)
subset_MACAU<-heatmap.2(as.matrix(subset_percents),na.color="grey",col=brewer.pal(9,"Reds"),trace=c("none"))

#make df
myDiffq_df<-data.frame(myDiffq)
head(myDiffq_df)
myDiffq_df$site<-paste(myDiffq_df$chr,myDiffq$start,myDiffq$start,sep=".")
head(myDiffq_df)

#myDiff20p
write.csv(myDiffq,file="RRBS_Sp_20samples_7of10_DMC_allCpG_qvalcutoff0.01_5x")

subset_methylkit<-heatmap.2(as.matrix(subset_percents),na.color="grey",col=brewer.pal(9,"Reds"),trace=c("none"))


###making a heatmap of the significant CpGs
setwd("~/clusteroutput/EpiMet_newgenome/filesforMACAU/")
sig_sites<-read.table("heritable50_sites_notrim60indv.txt", sep = '\t',header=F)
head(sig_sites)

#You can then do something like this
row.names(percents) %in% sig_sites
#This will give you a list of T/F values which you can then subset your table.  For example
subset_percents<-percents[row.names(percents) %in% sig_sites$V1, ]
subset_percents
nrow(subset_percents)
write.table(subset_percents,"MACAUoutput_notrim60indv_FDR_sigsites_percents", sep = "\t", row.names = F, col.names = T, quote=F)
####making a heatmap of methylkit sig sites

#You can then do something like this
row.names(percents) %in% myDiffq_df
#This will give you a list of T/F values which you can then subset your table.  For example

subset_percents<-percents[row.names(percents) %in% myDiffq_df$site, ]
head(subset_percents)
library(gplots)
library(RColorBrewer)
subset_methylkit<-heatmap.2(as.matrix(subset_percents),na.color="grey",col=brewer.pal(9,"Reds"),trace=c("none"))

#############LOOKING AT PERMUTATIONS##################################
####after MACAU is performed, read in the output table to do multiple testing correction
#source("https://bioconductor.org/biocLite.R")
#biocLite("qvalue")
library(qvalue)
setwd("~/clusteroutput/EpiMet_Sp_version2genome_released/MACAUfiles/permute/")
macau_out<-read.table("output_permute5_sp_allreads60indv.assoc.txt", sep = '\t',header=T)
head(macau_out)
macau_out$fdr <- lfdr(macau_out$pvalue)

write.table(macau_out,"MACAUoutput_permute5_sp_allreads60indv_FDR", sep = "\t", row.names = F, col.names = T, quote=F)


#########dendrogram#####################
library(ggplot2)
library(ggdendro)
library(gridExtra)
library(RColorBrewer)
library(colorRamps)
meth<-methylKit::unite(myobj,min.per.group=15L)
clus.obj<-clusterSamples(meth, dist = "manhattan", method="ward", plot=F)

# Assuming there is already a 'clus.obj' that contains the clustering
# of the samples

# Convert the clus.obj to a 'dendro' obj
dendro.data <- dendro_data(clus.obj)

# Build our data.frame for various annotation of the plots,
# Treatment, SampleIDs, Family, Family Shape and Color
df2 <- data.frame(
  treatment=factor(ifelse(grepl(pattern="S_", clus.obj$labels), "S", "T")),
  samples=factor(clus.obj$labels, levels=clus.obj$labels[clus.obj$order]),
  family=factor(unlist(strsplit(clus.obj$labels, split="_"))[seq(from=2, to=180, by=3)]),
  stringsAsFactors=F)

# This builds a set of combonations of shapes and colors for the family
# points.
shape.key <- data.frame(
  shape=rep(21:25, times=5), 
  colour=rep(c("A", "B", "C", "D", "E"), each=5), 
  stringsAsFactors=F)

df2$shape <- shape.key[as.integer(df2$family), 1]
df2$colour <- shape.key[as.integer(df2$family), 2]

# First Plot - Dendrogram with Family Markers and SampleIDs
# geom_segment to make the dendrogram
# geom_point to make the family markers
#
#   The scale_x_continous lets us adjust the labels, if you want to change
#   them, this would be spot, but you would need them in the right order.
#
#   The Colors for the family markers are set vial the scale_color_brewer,
#   If you want manual colors you can use scale_color_manual() instead.
#
#   For both plots, we're setting 'expand' in the x-axis to the same, this lines
#   the two plots together in the last step.
#
#   For geom_point, the positions are based on the number of rows in our df2
#   dataframe so if you want to use a different dataset, just need to define a
#   new df2 object where the number of rows is equal to the number of samples.
p1 <- ggplot(segment(dendro.data)) +
  scale_x_continuous(breaks=seq_along(dendro.data$labels$label),
                     labels=dendro.data$labels$label,
                     expand=c(0,0.6)) +
  scale_y_continuous() +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend)) +
  geom_point(data=df2[clus.obj$order, ], 
             aes(x=1:nrow(df2), y=rep(0, nrow(df2)), shape=shape, colour=colour, size=3, stroke=1)) +
  scale_color_brewer(palette="Dark2") +
  scale_shape_identity() +
  theme_dendro() +
  theme(axis.text.x=element_text(angle=90, hjust=1),
        legend.position="none",
        axis.text = element_text(size = 11))

# Second Plot - Plotting the Treatment Bar.
# geom_tile makes the actual bar
#
#   scale_fill_discrete allows us to adjust the labels of things although
#   you could also do this in the df2 data.frame
#
#   Adjusting the legend with both theme and guides. Theme adjusts position
#   while guide adjusts size of the colored rectangles in the legend.
#
#   Forcing spacing between the labels by adding extra spaces in the text.
p2 <- ggplot(df2, aes(samples, y=1, fill=factor(treatment))) + geom_tile() +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_discrete(expand=c(0,0.6)) +  
  scale_fill_discrete(name="Treatments",
                      labels=c("Stream  ", "Tank  ")) +
  theme(axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text=element_blank(),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        legend.position="bottom",
        legend.title=element_blank(),
        legend.text = element_text(size = 12)) +
  guides(fill=guide_legend(
    keywidth=0.1,
    keyheight=0.1,
    default.unit="inch"
  ))

# Finally putting the two plots together
#
# Not sure what a ggplotGrob is but its needed.
gp1 <- ggplotGrob(p1)
gp2 <- ggplotGrob(p2)
maxWidth = grid::unit.pmax(gp1$widths[2:5], gp2$widths[2:5])
gp1$widths[2:5] <- as.list(maxWidth)
gp2$widths[2:5] <- as.list(maxWidth)

# Set the layout and percent of space each plot will use and then draw the
# combined plot
grid.arrange(gp1, gp2, ncol=1, heights=c(7/8, 1/8))

plotdata<-ggplot_build(gp1)


###########MODIFY FOR USE##############
#HEATMAP - TISSUE
library(gplots)
####same issue with labeling as ABOVE
SpDMC<-read.table("MACAUfiles/output_EpiMet_sp_allreads/MACAUoutput_notrim60indv_FDR_sigsites_percents", header=T ,row.names=1, quote = NULL)

head(SpDMC)
#SpDMC$V1<-NULL

nrow(SpDMC)
#make variable for IDs of sample
tissue_heatmap<-c("T",
                  "T",
                  "T",
                  "S",
                  "S",
                  "S",
                  "T",
                  "T",
                  "T",
                  "T",
                  "S",
                  "S",
                  "S",
                  "T",
                  "T",
                  "S",
                  "S",
                  "S",
                  "T",
                  "T",
                  "S",
                  "S",
                  "T",
                  "S",
                  "T",
                  "T",
                  "T",
                  "T",
                  "S",
                  "T",
                  "T",
                  "S",
                  "S",
                  "S",
                  "T",
                  "S",
                  "T",
                  "T",
                  "T",
                  "T",
                  "S",
                  "S",
                  "T",
                  "S",
                  "S",
                  "S",
                  "S",
                  "S",
                  "S",
                  "S",
                  "T",
                  "S",
                  "S",
                  "S",
                  "S",
                  "T",
                  "T",
                  "T",
                  "T",
                  "T")
hclustfunc <- function(x) hclust(x, method="ward.")
col_anno <- as.matrix(ifelse(tissue_heatmap == "S", "#F8766D", "#00BFC4"))
colnames(col_anno) <- "Rearing"
CellDMR<-heatmap.2(as.matrix(SpDMC),na.color="grey",labRow=NA,labCol=NA,col=brewer.pal(9,"Reds"),trace=c("none"),ColSideColors=col_anno)

################getting total CpG covered and CpG covered at 10x for all samples########################



# Generating list of files
files <- list.files(pattern="*reformcov.txt")

# Building empty data.frame for storing data
data.coverage <- data.frame(Sample=character(length(files)),
                            Coverage=numeric(length(files)),
                            TotalCoverage=numeric(length(files)),
                            stringsAsFactors=F)

# Looping through the samples
for (i in 1:length(files)) {
  print(files[i])
  
  # Making Sample Name
  name <- paste0(gsub(pattern="_sp_allreads_reformcov.txt", replacement="", x=files[i]))    
  
  # Reading in the file
  tmp <- read.table(files[i], header=T)
  
  # Saving the data to the data.frame
  data.coverage$Sample[i] <- name
  data.coverage$Coverage[i] <- nrow(tmp[tmp$coverage >= 10, ])
  data.coverage$TotalCoverage[i] <- nrow(tmp)
}

write.table(data.coverage,"EpiMet_Sperm_metaCGcoverage.txt", sep = "\t", row.names = F, col.names = T, quote=F)




#closest gene to ALL CpG after uniting min 15 - did bedtools closest on the file, then imported and filtered to within 10k, then removed dupliates
#making bed of all combined
bed<-meth_df[,1:3]
head(bed)
write.table(bed,"EpiMet_Sp_allCpGunited15.bed", sep = "\t", row.names = F, col.names = F, quote=F)

#performed bedtools closest: ./software/bedtools2/bin/bedtools closest 
#-a /Users/mackenzie.gavery/clusteroutput/EpiMet_Sp_version2genome_released/EpiMet_Sp_allCpGunited15.sorted.bed 
#-b /Users/mackenzie.gavery/Desktop/genome_files/Omykiss_released_version2/bestblast/GeneLocation_50764gene_anddirection.sorted.bed 
#-D b > /Users/mackenzie.gavery/clusteroutput/EpiMet_Sp_version2genome_released/EpiMet_Sp_allCpGunited15.closestgene.txt 

#read in the output (some sites have 2 genes)
all_genes<-read.table("EpiMet_Sp_allCpGunited15.closestgene.txt", sep = '\t',header=F)
head(all_genes)

#filter to only sites with hits
library(dplyr)
all_genes_hitsonly<-filter(all_genes, V4 != ".")
head(all_genes_hitsonly)
nrow(all_genes_hitsonly)

#find absolute value of distance
all_genes_hitsonly$abs<-abs(all_genes_hitsonly$V10)
head(all_genes_hitsonly)
all_genes_win10k<-filter(all_genes_hitsonly, abs <= 10000)
nrow(all_genes_win10k)
head(all_genes_win10k)

#keep only unique genes
all_genes_win10k_nodupes<-all_genes_win10k %>% distinct(V7, .keep_all = TRUE)
nrow(all_genes_win10k_nodupes)
head(all_genes_win10k_nodupes)

#read in significant genes
sig_genes<-read.table("MACAUfiles/output_EpiMet_sp_allreads/EpiMet_Sp_DMCgenes_mockpvalue.txt", sep = '\t',header=T)
nrow(sig_genes)
#join ALL genes to sig genes
IPAtable<-left_join(all_genes_win10k_nodupes,sig_genes, by = c("V7" = "geneID"))
nrow(IPAtable)
head(IPAtable)
#new df with just gene and pvalue
IPAtable_minimize<-IPAtable[,c(7,12)]
head(IPAtable_minimize)
nrow(IPAtable_minimize)
#replace NA with 1
IPAtable_minimize[is.na(IPAtable_minimize)]<-"1"
#count how many zeros to make sure same number as sig genes
sum(IPAtable_minimize$pvalue == 0) 
#change column name
colnames(IPAtable_minimize)[1] <- "geneID"
#write the table
write.table(IPAtable_minimize,"EpiMet_Sp_DMC_IPAtable", sep = "\t", row.names = F, col.names = T, quote=F)
