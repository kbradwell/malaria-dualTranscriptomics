#source("https://bioconductor.org/biocLite.R")
#biocLite("fgsea")

library(ggplot2)
library(fgsea)
library(limma)
library(GEOquery)
library(Rcpp)
library(data.table)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packageVersion('fgsea') 
BiocManager::install("reactome.db")

source('https://bioconductor.org/biocLite.R')
biocLite('org.Hs.eg.db')
biocLite('org.Mm.eg.db')

library('org.Hs.eg.db')
library('org.Mm.eg.db')

columns(org.Hs.eg.db)

#gmt.file <- gmtPathways(system.file("extdata", "c5.bp.v7.1.symbols.gmt", package="fgsea"))
#pathways <- gmtPathways(gmt.file)
#str(head(pathways))

data(examplePathways)
data(exampleRanks)
exampleRanks
#fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=10000, maxSize=500)
# Testing only one pathway is implemented in a more efficient manner
fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=1000)
print(names(examplePathways))

setwd("/Users/kbradwell/Desktop")
getwd()

# # format the ranks for pheatmap (not run)
# rankdf <- read.csv(file = "./gsea_response_to_reviewers/human_TP_genesymbol_rankingstat_uniq.csv", header = TRUE, row.names = 1)
# rankdf
# ranksDE <- data.table(rankdf,keep.rownames = TRUE)
# ranksDE
# myranks <- ranksDE[order(RANK), list(rn, RANK)]
# myranks

GSEA = function(gene_list, GO_file) {
  ranks <- read.csv(file = gene_list, header = TRUE)
  #ranks <- read.table(ranks,header=TRUE, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$RANK, ranks$GENESYMBOL)
  str(ranks)
  
  # GO Pathways BP
  fgseaRes <- fgsea(GO_file, ranks, nperm=1000)
  
  #print(names(pathways.bp))
  # number of significant pathways at padj < 0.01
  sum(fgseaRes[, padj < 0.05])
  
  # plot the most significantly enriched pathway
  plotEnrichment(GO_file[[head(fgseaRes[order(pval), ], 1)$pathway]],ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway)
  
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(GO_file[topPathways], ranks, fgseaRes, gseaParam=0.5)
  
  output = list("Results" = fgseaRes)
  return(output)
}

pathways.bp <- gmtPathways("./gsea_response_to_reviewers/c5.bp.v7.1.symbols.gmt")
pathways.cc <- gmtPathways("./gsea_response_to_reviewers/c5.cc.v7.1.symbols.gmt")
pathways.mf <- gmtPathways("./gsea_response_to_reviewers/c5.mf.v7.1.symbols.gmt")
pathways.reactome <- gmtPathways("./gsea_response_to_reviewers/c2.cp.reactome.v7.1.symbols.gmt")

gene_list = "./gsea_response_to_reviewers/human_TP_genesymbol_rankingstat_uniq.csv"

res = GSEA(gene_list, pathways.bp)
dim(res$Results)
res$Plot

# reactome with reactomePathways()

GSEA_reactome = function(gene_list) {
  ranks <- read.csv(file = gene_list, header = TRUE)
  #ranks <- read.table(ranks,header=TRUE, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$RANK, ranks$GENESYMBOL)
  str(ranks)
  
  # convert gene symbols to Entrez gene IDs (for reactomePathways function)
  symbols <- c(names(ranks))
  symbols
  
  genemapping <- mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
  
  # format the new Entrez IDs to make sure any NAs or duplicated IDs are unique
  pseudoName <- 0
  entrezIDs <- c()
  for (i in c(1:length(genemapping))){
    pseudoName <- pseudoName + 1
    if(genemapping[[i]] %in% entrezIDs) {
      newID <- paste(genemapping[[i]], pseudoName, sep = "")
      entrezIDs <- c(entrezIDs,newID)
    } else {
      entrezIDs <- c(entrezIDs,genemapping[[i]])
    }
  }
  print(entrezIDs)
  print(length(entrezIDs))
  print(length(unique(entrezIDs)))
  
  ranks2 <- read.csv(file = gene_list, header = TRUE)
  ranks2$ENTREZ <- entrezIDs
  ranks2
  ranks2 <- subset( ranks2, select = -GENESYMBOL )
  ranks2
  ranks2 <- setNames(ranks2$RANK, ranks2$ENTREZ)
  str(ranks2)
  
  my_pathways <- reactomePathways(names(ranks2))
  
  # Reactome pathways have a median of 11 genes
  summary(sapply(my_pathways, length))
  
  fgsea_reactome <- fgsea(pathways = my_pathways, 
                          stats = ranks2,
                          minSize=15,
                          maxSize=500,
                          nperm=1000)
  
  head(fgsea_reactome[order(pval), ])
  
  sum(fgsea_reactome[, padj < 0.05])
  
  # plot the most significantly enriched pathway
  plotEnrichment(my_pathways[[head(fgsea_reactome[order(pval), ], 1)$pathway]],ranks2) + labs(title=head(fgsea_reactome[order(pval), ], 1)$pathway)
  
  topPathwaysUp <- fgsea_reactome[ES > 0][head(order(pval), n=15), pathway]
  topPathwaysDown <- fgsea_reactome[ES < 0][head(order(pval), n=15), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  #topPathways
  plotGseaTable(my_pathways[topPathways], ranks2, fgsea_reactome, gseaParam=0.5,colwidths = c(7.2, 1.8, 0.6, 0.7, 0.7))
  
  output = list("Results" = fgsea_reactome)
  return(output)
}

#gene_list = "./gsea_response_to_reviewers/human_TP_genesymbol_rankingstat_uniq.csv"
gene_list = "./gsea_response_to_reviewers/human_P1P2_genesymbol_rankingstat_uniq.csv"
#gene_list = "./gsea_response_to_reviewers/human_P1P3_genesymbol_rankingstat_uniq.csv"

set.seed(42)
res = GSEA_reactome(gene_list)
dim(res$Results)
res$Results
print(sum(res$Results[, padj < 0.05]))
#fwrite(res$Results, file="fgseaRes_reactome_human_TP.txt", sep="\t", sep2=c("", " ", ""))
fwrite(res$Results, file="fgseaRes_reactome_human_P1P2.txt", sep="\t", sep2=c("", " ", ""))

# # can't fund function mapIdsList
# res$Results[, leadingEdge := mapIdsList(
#   x=org.Hs.eg.db, 
#   keys=leadingEdge,
#   keytype="ENTREZID", 
#   column="SYMBOL")]
# fwrite(res$Results, file="fgseaRes_reactome_human_TP_leading_edge.txt", sep="\t", sep2=c("", " ", ""))