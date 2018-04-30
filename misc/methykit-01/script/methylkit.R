
library(methylKit)

file.list_a=list('data/zr2096_1_dedup.sorted.bam',
 'data/zr2096_2_dedup.sorted.bam',
  'data/zr2096_3_dedup.sorted.bam',
  'data/zr2096_4_dedup.sorted.bam',
  'data/zr2096_5_dedup.sorted.bam',
  'data/zr2096_6_dedup.sorted.bam',
  'data/zr2096_7_dedup.sorted.bam',
  'data/zr2096_8_dedup.sorted.bam',
  'data/zr2096_9_dedup.sorted.bam',
  'data/zr2096_10_dedup.sorted.bam'
  )



myobj2 = processBismarkAln(location = file.list_a, sample.id = list("1","2","3","4","5","6","7","8","9","10"), assembly = "v3", read.context="CpG", mincov=10, treatment = c(0,0,0,0,0,1,1,1,1,1))



getMethylationStats(myobj2[[1]],plot=FALSE,both.strands=FALSE)

getMethylationStats(myobj2[[1]],plot=TRUE,both.strands=FALSE)

getCoverageStats(myobj2[[10]],plot=TRUE,both.strands=FALSE)




meth=unite(myobj2)

head(meth)



getCorrelation(meth,plot=TRUE)


clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)



hc = clusterSamples(meth, dist="correlation", method="ward", plot=FALSE)

PCASamples(meth, screeplot=TRUE)




PCASamples(meth)


myDiff=calculateDiffMeth(meth)



myDiff25p=getMethylDiff(myDiff,difference=25,qvalue=0.01)


myDiff50p <- getMethylDiff(myDiff,difference=50,qvalue=0.01)

write.table(myDiff50p, file = "analyses/myDiff50p.tab")

View(myDiff50p)


write.csv(myDiff25p, file = "analyses/myDiff25p.csv")

