library(tidyverse)

file.name <- "5pCtoT_freq_2B_01.txt"
df<- read.delim(file.name)

head(df)

ggplot()+
  geom_line(data = df, aes(x= pos, y = X5pC.T), color = "red")+
  ylim(0,0.30)+
  ylab('Frequency of C to T')+
  xlab("distance from 5' end (bp)")+
  theme_bw()

