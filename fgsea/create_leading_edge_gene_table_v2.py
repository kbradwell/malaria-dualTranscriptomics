from sys import argv

script, mappingfile, resfile, outfile, comparison = argv

MAP_IFH = open(mappingfile,'r') # ncbiID2ncbiDesc.txt
GSEARES_IFH = open(resfile,'r') # GSEA results file
OFH = open(outfile,'w') # table of leading edge genes per pathway

comparisonStr = str(comparison)

mappingHash = {}
lines = MAP_IFH.readlines()[1:]
for l in lines:
	cols = l.strip().split("\t")
	if len(cols)>1:
		ncbiID = cols[0]
		ncbiDesc = cols[1] # formerly known as Entrez gene ID
		if ncbiID not in mappingHash:
			mappingHash[ncbiID] = ncbiDesc
		else:
			print("duplicated")
			print(ncbiID)

#for entrez in mappingHash:
#	print(mappingHash[entrez])

OFH.write("comparison\tpathway\tpval\tpadj\tES\tNES\tnMoreExtreme\tsize\tleading edge gene NCBI ID\tgene description\n")
lines = GSEARES_IFH.readlines()[1:]

NESlines = []
for l in lines:
	cols = l.strip().split("\t")
	nes = float(cols[4])
	linetup = (nes,l)
	NESlines.append(linetup)

sortedNES = sorted(NESlines, key = lambda x: x[0], reverse=True)
sortedLines = [x[1] for x in sortedNES]

for l in sortedLines:
	cols = l.strip().split("\t")
	padj = float(cols[2])
	leadingEdge = cols[7]
	leadingEdgeList = leadingEdge.strip().split(" ")
	OFH.write(comparisonStr)
	OFH.write("\t")
	enrichmentInfo = "\t".join(cols[0:7])
	OFH.write(enrichmentInfo)
	if padj < 0.05:
		geneCount = 0
		for gene in leadingEdgeList:
			if geneCount == 0:
				OFH.write("\t")
				OFH.write(gene)
				OFH.write("\t")
				OFH.write(mappingHash[gene])
				OFH.write("\n")
				geneCount = 1
			else:
				OFH.write("\t\t\t\t\t\t\t\t")
				OFH.write(gene)
				OFH.write("\t")
				OFH.write(mappingHash[gene])
				OFH.write("\n")
	else:
		OFH.write("\t--\t--\n")
	
MAP_IFH.close()
GSEARES_IFH.close()
OFH.close()

	
