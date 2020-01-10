# differential expression for Mali Paxgene samples

# source("http://bioconductor.org/biocLite.R")
# biocLite(pkgs=c("Rsubread","limma","edgeR"))

library(Rsubread)
library(edgeR)
library(ggplot2) 
library(ggrepel)
library(gridExtra)

#### INPUT FILES ##########################################################################

# readsDF <- read.table("./str_pysamcounts_v6_hg3_noglobin_autosomes_libsizenorm.txt",header = TRUE,sep="\t")
readsDF <- read.table("./str_pysamcountsv6_pf.libsizenorm.txt",header = TRUE,sep="\t")

###########################################################################################

dfcounts <- as.data.frame(readsDF)
annotLocs <- grep("GeneID|Length", colnames(dfcounts))
grp1Locs <- grep("_patient1", colnames(dfcounts))
grp2Locs <- grep("_patient2", colnames(dfcounts))
grp3Locs <- grep("_patient3", colnames(dfcounts))
x <- dfcounts[,c(annotLocs, grp1Locs, grp2Locs, grp3Locs)]

colnames(x)
head(x, 15)

# drop the "Length" column from the data frame - the only columns present should just contain counts
row.names(x) <- x$GeneID
cat("Row names assigned\n")
drops <- c("GeneID","Length")
x <- x[ , !(names(x) %in% drops)]

colnames(x)
head(x, 15)

group <- factor(c(1,1,1,1,1,2,2,2,2,2,3,3,3,3,3))
timepoints <- c(1,2,3,4,5,1,2,3,4,5,1,2,3,4,5)
y <- DGEList(counts=x)

# filter out non-informative genes (low counts) - must have >=X cpm in Y or more samples
keep <- rowSums(cpm(y)>10) >= 7
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]


y <- calcNormFactors(y)

# add in stage info for Pf
stage_estimates <- read.table("./pf_stages2.txt",header = TRUE,sep="\t", row.names = 1)
CIBER_Ring_stage <- stage_estimates$Ring
CIBER_Troph_stage <- stage_estimates$Troph30
# design <- model.matrix(~group+timepoints)
design <- model.matrix(~group+timepoints+CIBER_Ring_stage)
# rownames(design) <- colnames(y)
# DEgenes <- estimateDisp(y,design)

#### Tagwise dispersion ####

cat("Tagwise!\n")
DEgenes <- estimateGLMCommonDisp(y, design)
# To estimate trended dispersions:
# y <- estimateGLMTrendedDisp(DEgenes, design)
# To estimate tagwise dispersions:
DEgenes <- estimateGLMTagwiseDisp(DEgenes, design)

############################

cat("Design columns\n")
colnames(design)
cat("Design\n")
design

# DEgenes$common.dispersion

cat("Dispersion estimated\n")

# perform quasi-likelihood F-tests
fit <- glmQLFit(DEgenes,design)
cat("fit columns\n")
colnames(fit)

qlf <- glmQLFTest(fit,coef=4)
# qlf <- glmQLFTest(fit,coef=2:3)
cat("The total number of genes significantly up-regulated or down-regulated at 5% FDR\n")
summary(decideTests(qlf$table$PValue))
#topTags(qlf,n=100)
num2show = nrow( qlf$table )
topTags(qlf,n = num2show )
	
head(sort(round(p.adjust(qlf$table$PValue, "BH"), 3)),5)


