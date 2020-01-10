import sys
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl

'''
IQR and Index of Dispersion - based filtering of read counts
'''

script, counts, outfile = argv

cutoffIoD = 0.1 # cutoff for IoD
cutoffIQR = 0.3 # cutoff for IQR

IFH=open(counts,'r')
OFH=open(outfile,'w')

tooLow=0
goodDispersion=0
firstLine=0
for line in IFH.readlines():
	if firstLine==0:
		OFH.write(line)
		firstLine+=1
	else:
		cols=line.split("\t")
		counts = cols[2:]
		floatCounts = ([float(xi) for xi in counts]) 
		# calculate variance
		resIoD=(np.var(floatCounts)/np.mean(floatCounts))
		resIQR= np.subtract(*np.percentile(floatCounts, [75, 25]))
		if (resIoD>cutoffIoD and resIQR>cutoffIQR):
			OFH.write(line)
			goodDispersion+=1
		else:
			tooLow+=1
			print line

print "Number of genes below threshold:", tooLow
print "Number of genes above threshold (kept):", goodDispersion

IFH.close()
OFH.close()
