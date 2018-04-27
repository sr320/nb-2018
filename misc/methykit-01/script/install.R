install.packages("devtools")

library(devtools)

source("https://bioconductor.org/biocLite.R")
biocLite("methylKit")

install_github("al2na/methylKit", build_vignettes=FALSE, 
               repos=BiocInstaller::biocinstallRepos(),
               dependencies=TRUE)


library(methylKit)

