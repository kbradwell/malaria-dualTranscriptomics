import glob, os, re
os.chdir("./")

import sys
from sys import argv
import re

'''
Get reference allele frequencies from pileup files
'''

script, excluderegions, outfile = argv

IFH=open(excluderegions,'r')
OFH=open(outfile,'w')

import glob, os, re
os.chdir("./")

##########################################################################################
#
# SAVE EXCLUSION SITE INFO
#
##########################################################################################

exclusionList={}
for line in IFH:
	cols=line.split('\t')
	chr = cols[0]
	start = cols[3]
	end = cols[4]
	for x in range(int(start),int(end)+1):
		if chr not in exclusionList:
			exclusionList[chr]={}
		exclusionList[chr][str(x)]=1

IFH.close()

##########################################################################################
#
# PRINT RAF PER SITE
#
##########################################################################################

skipStrings = (">", "<")
for file in glob.glob("*.bam.pileup"):
    f = open(file, "r")
    for line in f:
    	cols=line.split('\t')
    	bases = cols[4]
    	refAllele = cols[2]
    	cov = int(cols[3])
    	chr = cols[0]
    	pos = cols[1]
    	loc = "_".join([chr,pos])
    	#print bases
    	if any(s in bases for s in skipStrings):
    		pass
    	elif "N" in refAllele:
    		pass
    	elif cov <50:
    		pass
    	elif ((chr in exclusionList) and (str(pos) in exclusionList[chr])):
    		pass
    	else:
    		RA = len(re.findall("\.|,", bases))
    		RAF = RA/float(cov)
    		OFH.write(loc)
    		OFH.write("\t")
    		OFH.write("{0:.2f}".format(RAF))
    		OFH.write("\t")
    		OFH.write(str(file))
    		OFH.write("\n")
OFH.close()
