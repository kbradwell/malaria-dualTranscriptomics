from sys import argv
import re

script, globinlist, countsIN, countsOUT = argv

globins=open(globinlist,'r')
IFH=open(countsIN,'r')
OFH=open(countsOUT,'w')

globinGenes=[]
for line in globins:
	gene=(line.strip()).replace("\"","")
	globinGenes.append(gene)

#print globinGenes

for line in IFH:
	genename=line.split("\t")[0]
	if genename.strip() not in globinGenes:
		OFH.write(line)
	