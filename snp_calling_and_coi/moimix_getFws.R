# have to run this interactively on grid (qlogin -> cd dir -> use r-3.3.1 -> R -> copy/paste blocks of commands)

library(SeqArray)
library(moimix)
seqVCF2GDS("./pf215_gatk_output_split_filt_maxallele2_maxmissing0.8_snps-only.vcf.recode.vcf.gz", "./pf_gatk_split.gds")

isolates <- seqOpen("./pf_gatk_split.gds")
seqSummary(isolates)

# save sample identifiers
sample.id <- seqGetData(isolates, "sample.id")

# get genomic coordinates of all variants
#coords <- getCoordinates(isolates)
#head(coords)

fws_all <- getFws(isolates)

#hist(fws_all)

# see if our sample that we estimated is multiclonal according to fws
#fws_all["./uniq_bothR1R2mapped_primary_hisat2_pf_36D.namesort.fixmate.sort.nodup.bam.pileup"] < 0.95

for (samp in fws_all){
	print(paste("Fws", samp))
}

for (samp in sample.id){
	print(paste("Fws", samp, fws_all[samp]))
}
