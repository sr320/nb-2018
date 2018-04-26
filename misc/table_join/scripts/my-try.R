

#importing tables - I edited column title in as odd character present
qubit <- read.csv("data/Qubit-consolidated.csv", sep=",")


RNA <- read.csv("data/samples-for-RNA-isolation.csv", sep=",")


#loading dplyr package 
library(dplyr)

#inner join
crab_data <- inner_join(qubit, RNA, by = "tube_number")

#view
crab_data
