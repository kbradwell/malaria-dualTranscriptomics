library(ggplot2)
#Stacked + percent
ggplot(multispp_all_counts, aes(fill=SPP, y=NUMREADS, x=SAMP)) + geom_bar(position="fill", stat="identity") + xlab("sample") + ylab("Proportion of reads") + labs(fill = "species")