1. permutations_mali_v2 obtains the permutations on the original read counts
2. split_csv is then used on the original and permuted datasets to allow faster (parallel) computation of read count correlations downstream
3. readcounts2spearman_v4 calculates the correlations on the split original and permuted readcounts
4. get_permutation_pvals_v2 allows for p-value/FDR-type calculations on the correlation results
