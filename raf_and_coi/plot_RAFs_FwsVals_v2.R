library(ggplot2)
library(grDevices)

MALIRAFs50Xv2 <- read.csv(file="./outRAFs50X_noMGFs_v2_patientIDs.csv", header=TRUE, sep=",")

pdf("RAF_hist_FwsVals_patientIDs.pdf", 12, 9)

sampVals <- c("1A", "1B","1C","1D","1E","2A","2B","2C","2D","2E","3A","3B","3C","3D","3E")
fwsVals <- c("0.86", "0.98", "0.80", "0.98", "0.96", "0.99", "0.99", "0.73", "0.99", "0.99", "0.82", "0.83", "0.99", "0.85", "0.98")


# create a list of labels using bquotw
# labs <- Map(.beta = sampVals, .p = fwsVals, f = function(.beta,.p) bquote(list(beta == .(.beta), italic(p) == .(.p))))
labs <- Map(.p = fwsVals, f = function(.p) bquote(list(italic(F[WS]) == .(.p))))
# coerce to a character representation for parse=TRUE to work within 
len <-length(levels(MALIRAFs50Xv2$SAMP))
vars <- data.frame(expand.grid(levels(MALIRAFs50Xv2$SAMP)))
colnames(vars) <- c("SAMP")

# geom_text
dat <- data.frame(x = rep(0.5, len), y = rep(150, len), vars, labels = sapply(labs,deparse))


p <- ggplot(data = MALIRAFs50Xv2, aes(x=RAF)) + geom_histogram(color="black", fill="white", alpha = 0.4,bins=80) + theme_bw() + theme(panel.spacing = unit(1.5, "lines"),axis.text = element_text(size = 12),axis.title.x = element_text(size = 16),axis.title.y = element_text(size = 16),strip.text.x = element_text(size = 14)) + coord_cartesian(ylim=c(0,200)) + ylab("count (number of loci)") + xlab("Reference Allele Frequency") #ylim(c(0, 500))

histplot <- p + facet_wrap( ~ SAMP, ncol = 5) + geom_text(data = dat, size = 6, aes(x=x,y=y,label=labels), parse=TRUE)

plot(histplot)

dev.off()