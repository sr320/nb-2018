

library(dplyr)

#read in data
crab_big3 <- read.csv("https://raw.githubusercontent.com/RobertsLab/project-crab/master/data/goodsamples.csv")



#It contains all the good samples (20 ng/µL or more RNA). I am wondering if there is a way to sort the rows based on treatment (warm/cold/ambient) and infection status (0/ 1) so that I can more easily visualize what I have within each category.


group_by(crab_big3, temperature_treatment) %>%
  summarize(count = n()) 

#output
# A tibble: 3 x 2
#temperature_treatment count
#<fct>                 <int>
#1 ambient                  37
#2 cold                     33
#3 warm                      3

# filter data for ambient and inf = 1 and counts
filter(crab_big3, temperature_treatment == "ambient", infection_status == 1) %>%
  summarize(count = n()) 

# filter data for ambient and inf = 0 and counts
filter(crab_big3, temperature_treatment == "ambient", infection_status == 0) %>%
  summarize(count = n()) 

library(ggplot2)



#not pretty but histogram for easy look
ggplot(crab_big3, aes(x = Sample_Day)) +
  geom_bar(aes(fill = temperature_treatment)) +
  facet_wrap(~infection_status)

#not pretty but histogram for easy look
ggplot(crab_big3, aes(x = infection_status)) +
  geom_bar(aes(fill = temperature_treatment)) +
  facet_wrap(~Sample_Day)+
  scale_y_continuous(breaks = round(seq(min(0), max(15), by = 1),1))



