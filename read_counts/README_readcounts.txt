1. use add2gtf_ensemblID to add ENSEMBL IDs to gtf if they aren't present already (read count will be performed for each unique genomic loci via ENSEMBL ID)
2. pysamcounts_mp_stranded_unstr_v6 performs readcounts using multiprocessing in Python
3. remove globin count rows using remove_globin_count_rows
4. normalize counts with apply_libsize_scalingfactor (command is given in run_normalize_count)
