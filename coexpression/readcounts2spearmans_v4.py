import sys
from sys import argv
import glob
import os
#import rpy2
#sys.path.append("/usr/local/packages/Python-2.7/lib/python2.7/site-packages/scipy/stats")
#print sys.path
from scipy.stats import spearmanr
#import scipy

script, file_name = argv

pffile = open('./split_pf.csv','r')

infilename = os.path.basename(file_name)
#print(file_name)
hgCounts={}
pfCounts={}
hgfile = open(file_name, 'r')

for line in hgfile:
	cols=line.strip().split(",")
	genename=cols[0]
	hgCounts[genename]=[]
	hgCounts[genename]=[float(i) for i in cols[1:]]

for line in pffile:
	cols=line.strip().split(",")
	genename=cols[0]
	pfCounts[genename]=[]
	pfCounts[genename]=[float(i) for i in cols[1:]]

belowThresh=0
outfile="correlations_"+str(infilename)+".out"
OFH=open(outfile,'w')
outfile2="low_correlations_"+str(infilename)+".out"
OFH2=open(outfile2,'w')
for humangene in hgCounts:
	for parasitegene in pfCounts:
		uniqID=humangene+"_"+parasitegene
		corr, p_value = spearmanr(hgCounts[humangene],pfCounts[parasitegene])
		if abs(corr)>0.8:
			OFH.write(str(infilename))
			OFH.write("\t")
			OFH.write(uniqID)
			OFH.write("\t")
			OFH.write(str(corr))
			OFH.write("\t")
			OFH.write(str(p_value))
			OFH.write("\n")
		else:
			belowThresh+=1
lessSig = str(infilename)+","+"<"+","+str(belowThresh)
OFH2.write(lessSig)
OFH2.write("\n")

hgfile.close()
pffile.close()
OFH.close()
OFH2.close()
