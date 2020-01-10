#library(tidyverse)

hg_pf_TPM_filtIQRIoD <- read.csv("./hg_pf_TPM_filtIQRIoD_out.csv")

dim(hg_pf_TPM_filtIQRIoD)
head(hg_pf_TPM_filtIQRIoD)

df <- hg_pf_TPM_filtIQRIoD
head(df)
# get the hg and pf subsets of the dataframe
hgDF <- df[c(1:12487),]
pfDF <- df[c(12488:17397),]

head(hgDF)
tail(hgDF)
head(pfDF)
tail(pfDF)

sampCols <- 2:16

for(p in 1:500) {
  colnames(hgDF) = colnames(df)
  
  pfDF2 <- pfDF[,c(1:1,sample(sampCols,size = 15, replace = FALSE))]
  print(head(pfDF2))
  colnames(pfDF2) = colnames(df)
  
  # merge the datasets together and use this new df for all remaining processing
  hgpfDF2 <- rbind(hgDF,pfDF2)
  #print(tail(hgpfDF2))
  outfile <- paste(c("expr_",p,"_perm_out.csv"), collapse = "")
  print(outfile)
  write.table(hgpfDF2, outfile, sep=",",  quote=FALSE, row.names=F)

} # end of loop for randomization