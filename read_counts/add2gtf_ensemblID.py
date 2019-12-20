from sys import argv
import re

# add on the Ensembl IDs in gtf files

script, lookupfile, gtffile, outgtf = argv

IFH = open(lookupfile,'r')
IFH2 = open(gtffile,'r')
OFH = open(outgtf,'w')

ensemblHash={}
# read the biomart info on NM_ to EnsemblID conversion (ensembl_ids_NM-ONLY.txt)
for line in IFH:
	ensemblID=line.split("\t")[0]
	refSeqNM_ID=(line.split("\t")[1]).strip()
	if refSeqNM_ID not in ensemblHash:
		ensemblHash[refSeqNM_ID]=[] # multiple mRNAs can map to on gene, so many values may be the same
		ensemblHash[refSeqNM_ID].append(ensemblID)
	else:
		if ensemblID not in ensemblHash[refSeqNM_ID]:
			ensemblHash[refSeqNM_ID].append(ensemblID)

notFoundCount=0
foundCount=0
multiCount=0
# add the EnsemblID information into the GTF file
for line in IFH2:
	if not line.startswith("##"):
		if "exon" in line:
			toKeepPart1 = (line.split(";")[0]).replace("gene_id","transcript_id") + ";"
			toKeepPart2 = (line.split(";")[1]).replace("geneName","gene_name") + ";"
			# add the ensembl_id as the middle field for the last column of the GTF file
			NMnum=(((toKeepPart1.split("\"")[1]).replace("\"","")).replace(";","")).split(".")[0].strip()
			#print NMnum
			if NMnum in ensemblHash:
				if len(ensemblHash[NMnum])==1:
					ensembl="\"" + str(ensemblHash[NMnum][0]) + "\"" + ";"
					newgtfLine = toKeepPart1 + " gene_id " + ensembl + toKeepPart2
					OFH.write(newgtfLine)
					OFH.write("\n")
					foundCount+=1
				else:
					multiCount+=1
			else:
				notFoundCount+=1
			
		else:
			OFH.write(line)
			
print notFoundCount
print foundCount
print multiCount
	
IFH.close()
IFH2.close()
OFH.close()
