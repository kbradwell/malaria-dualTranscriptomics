# Malaria dual transcriptomics

Custom scripts for the analysis of host-parasite dual RNA-Seq data. 

Directories included in this repository:  

## check_for_multispp

Map reads to different Plasmodium genomes concatenated in a single fasta file and plot percentage of reads to each genome

## read_counts

Obtaining read counts per gene and normalization of counts over samples

## differential_expression

edgeR script to determine differential expression over individuals and number of infections

## fgsea

Run pre-ranked GSEA on human and parasite genes

## coexpression

Spearmans correlations and permutations to determine whether host genes and parasite genes are coexpressed at a greater rate than expected by chance

## raf_and_coi

Analysis of reference allele frequencies and complexity of infection
