

#loading library
library(methylKit)


#copying dedup sorted bams from Bismark to data dir
#output on macbook: /Users/sr320/Dropbox/wd/18-04-27
#these are from subset (not full data set)
#creating list of files
file.list=list(  'data/zr2096_1_dedup.sorted.bam',
                 'data/zr2096_2_dedup.sorted.bam',
                 'data/zr2096_3_dedup.sorted.bam',
                 'data/zr2096_4_dedup.sorted.bam',
                 'data/zr2096_5_dedup.sorted.bam',
                 'data/zr2096_6_dedup.sorted.bam',
                 'data/zr2096_6_dedup.sorted.bam',
                 'data/zr2096_8_dedup.sorted.bam',
                 'data/zr2096_9_dedup.sorted.bam',
                 'data/zr2096_10_dedup.sorted.bam'
)


#bringing in bam files , setting context, AND coverage
myobj = processBismarkAln(location = file.list, sample.id = list("1","2","3","4","5","6","7","8","9","10"), assembly = "v3", read.context="CpG", mincov=3, treatment = c(0,0,0,0,0,1,1,1,1,1))



getMethylationStats(myobj[[1]],plot=FALSE,both.strands=FALSE)

getMethylationStats(myobj[[1]],plot=FALSE,both.strands=TRUE)



getMethylationStats(myobj[[4]],plot=TRUE,both.strands=FALSE)

getCoverageStats(myobj[[4]],plot=TRUE,both.strands=FALSE)




meth=unite(myobj)



getCorrelation(meth,plot=TRUE)


clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)




hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)

PCASamples(meth, screeplot=TRUE)



#PCA
PCASamples(meth)

PCASamples(meth,screeplot=FALSE, adj.lim=c(0.0004,0.1),
           scale=TRUE,center=TRUE,comp=c(1,2),transpose=TRUE,sd.filter=TRUE,
           sd.threshold=0.5,filterByQuantile=TRUE,obj.return=FALSE)


#DML gene list
myDiff=calculateDiffMeth(meth)



myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)


myDiff50p=getMethylDiff(myDiff,difference=50,qvalue=0.01)


# get hyper-methylated
hyper=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hyper")

hyper

# get hypo-methylated
hypo=getMethylDiff(myDiff,difference=25,qvalue=0.01,type="hypo")



# Examples


data(methylKit)

# get differentially methylated bases/regions with specific cutoffs
all.diff=getMethylDiff(methylDiff.obj,difference=25,qvalue=0.01,type="all")

# get hyper-methylated
hyper=getMethylDiff(methylDiff.obj,difference=25,qvalue=0.01,type="hyper")

# get hypo-methylated
hypo=getMethylDiff(methylDiff.obj,difference=25,qvalue=0.01,type="hypo")