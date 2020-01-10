from sys import argv
import glob
import os
import numpy

script, file_name, outfile = argv

list_of_files = glob.glob('./expr*out.csv_correlations.txt')
print len(list_of_files)

r2hash = {}

for i in numpy.arange(0.8,1,0.01):
	r2hash[i]={}
	r2hash[i]["permlist"]=[]
	r2hash[i]["obs"]=0

#print r2hash

# save the observed number of pairs for each R2 threshold
print "saving observed pairs counts"
observations = open(file_name,'r')
for line in observations:
	cols = line. split()
	r2 = float(cols[2])
	for threshold in r2hash:
		if abs(r2)>=threshold:
			r2hash[threshold]["obs"]+=1

observations.close()

print "saving empirical pairs counts"
for f in list_of_files:
	PERM = open(f,'r')
	numPairs = {}
	print str(f)
	for i in numpy.arange(0.8,1,0.01):
		numPairs[i]=0
	for line in PERM:
		cols = line. split()
		r2 = float(cols[2])
		for threshold in r2hash:
			if abs(r2)>=threshold:
				numPairs[threshold]+=1
	for t in numPairs:
		r2hash[t]["permlist"].append(numPairs[t])
	PERM.close()
	
OFH = open(outfile,'w')

print "saving p-values"
OFH.write("R2\tObs#Pairs\tMean#PermutationPairs\t")
for i in range(1,len(list_of_files)+1):
	permutation="PERM"+str(i)+"\t"
	OFH.write(permutation)
OFH.write("total>=Obs#Pairs\tpval\n")
for r2thresh in sorted(r2hash):
	avgPermPairs=sum(r2hash[r2thresh]["permlist"])/float(len(r2hash[r2thresh]["permlist"]))
	r2hash[r2thresh]["permMeanPairs"]=avgPermPairs
	aboveObsPairs=sum(i >= r2hash[r2thresh]["obs"] for i in r2hash[r2thresh]["permlist"])
	pval=aboveObsPairs/float(len(r2hash[r2thresh]["permlist"]))
	OFH.write(str(r2thresh))
	OFH.write("\t")
	OFH.write(str(r2hash[r2thresh]["obs"]))
	OFH.write("\t")
	OFH.write(str(avgPermPairs))
	OFH.write("\t")
	for i in r2hash[r2thresh]["permlist"]:
		OFH.write(str(i))
		OFH.write("\t")	
	OFH.write(str(aboveObsPairs))
	OFH.write("\t")
	OFH.write(str(pval))
	OFH.write("\n")

OFH.close()

				
