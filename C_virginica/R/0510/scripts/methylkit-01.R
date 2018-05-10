

#loading library
library(methylKit)


#copying dedup sorted bams from Bismark to data dir
#output on macbook: /Users/sr320/Dropbox/wd/18-04-27

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

getMethylationStats(myobj2[[1]],plot=TRUE,both.strands=FALSE)

getCoverageStats(myobj2[[4]],plot=TRUE,both.strands=FALSE)




meth=unite(myobj2)



getCorrelation(meth,plot=TRUE)


clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)