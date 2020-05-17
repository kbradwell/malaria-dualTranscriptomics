human: 
1. run gsea_human.r (I used the second function in this script which uses reactomePathways(), but if you want to use ENSEMBL ID directly then use the first function) 
2. run create_leading_edge_gene_table_v2.py to produce a formatted results table (comparison is e.g. "infection_number" or "patientX_v_patientY" and will be converted to a string from the command line
parasite: 
1. create a Plasmodium-specific gmt file with create_parasite_gmt_file.py
2. run gsea_parasite.r
