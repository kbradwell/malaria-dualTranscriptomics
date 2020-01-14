
library(SeqArray)
library(moimix)
seqVCF2GDS("./pf.vcf.recode.vcf.gz", "./pf_split.gds")

isolates <- seqOpen("./pf_split.gds")
seqSummary(isolates)

# save sample identifiers
sample.id <- seqGetData(isolates, "sample.id")

# get genomic coordinates of all variants
#coords <- getCoordinates(isolates)
#head(coords)

fws_all <- getFws(isolates)

#hist(fws_all)

for (samp in fws_all){
	print(paste("Fws", samp))
}

for (samp in sample.id){
	print(paste("Fws", samp, fws_all[samp]))
}
