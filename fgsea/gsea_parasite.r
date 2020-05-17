#remotes::update_packages("rlang")
install.packages("dplyr")
library(ggplot2)
library(fgsea)
library(limma)
library(GEOquery)
library(Rcpp)
library(data.table)
library('org.Hs.eg.db')

#install.packages("msigdbr")
#install.packages("Rcpp")

#library(msigdbr)

setwd("/Users/kbradwell/Desktop")
getwd()

# msigdbr_show_species() doesn't have Plasmodium

GSEA = function(gene_list, GO_file) {
  ranks <- read.csv(file = gene_list, header = TRUE)
  ranks
  #ranks <- read.table(ranks,header=TRUE, colClasses = c("character", "numeric"))
  ranks <- setNames(ranks$RANK, ranks$GENEID)
  str(ranks)
  
  # Reactome
  fgseaRes <- fgsea(GO_file, ranks, nperm=1000)
  
  #print(names(pathways.reactome))
  # number of significant pathways at padj < 0.01
  sum(fgseaRes[, padj < 0.5])
  
  # plot the most significantly enriched pathway
  plotEnrichment(GO_file[[head(fgseaRes[order(pval), ], 1)$pathway]],ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway)
  
  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]
  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  plotGseaTable(GO_file[topPathways], ranks, fgseaRes, gseaParam=0.5)
  
  output = list("Results" = fgseaRes)
  return(output)
  
}

pathways.reactome <- gmtPathways("./gsea_response_to_reviewers/plasmodium.gmt")

#gene_list = "./gsea_response_to_reviewers/plasmodium_TP_geneID_rankingstat.csv"
#gene_list = "./gsea_response_to_reviewers/plasmodium_P1P2_geneID_rankingstat.csv"
gene_list = "./gsea_response_to_reviewers/plasmodium_P1P3_geneID_rankingstat.csv"

res = GSEA(gene_list, pathways.reactome)
dim(res$Results)
print(sum(res$Results[, padj < 0.3]))


