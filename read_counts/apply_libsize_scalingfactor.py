from sys import argv
import re

# normalize gene counts to match the lib size of the smallest sample

script, countsIN, countsOUT = argv

OFH=open(countsOUT,'w')

# save the lib size (sum of all read counts) for each sample
with open(countsIN) as f:
	first_line = f.readline()
	OFH.write(first_line)
	samplesList = first_line.strip().split("\t")
	samplesIndex=0
	samplesIndexLookup={}
	sampleLibCount={}
	linesList=[]
	for s in samplesList[2:]:
		samp=s.strip()+"_"+str(samplesList.index(s))
		samplesIndexLookup[samp]=samplesIndex
		samplesIndex+=1
	for line in f.readlines():
		linesList.append(line.strip())
		counts=line.strip().split("\t")[2:]
		for s in samplesIndexLookup:
			if s not in sampleLibCount:
				sampleLibCount[s]=0
				sampleLibCount[s]+=int(float(counts[samplesIndexLookup[s]]))
			else:
				sampleLibCount[s]+=int(float(counts[samplesIndexLookup[s]]))

minLibSize=min(sampleLibCount[s] for s in sampleLibCount)
#print minLibSize

#calculate scaling factor required for each sample
scalingFactors={}
for s in sorted(sampleLibCount):
	fac=minLibSize/float(sampleLibCount[s])
	print s, fac, sampleLibCount[s]
	scalingFactors[s]=fac

#write new counts to file, applying the scaling factors
newLibSizes={}
for l in linesList:
	#OFH.write(l)
	desc=l.split("\t")[0:2]
	desc2print="\t".join(desc)
	#print desc2print
	OFH.write(desc2print)
	OFH.write("\t")
	counts=l.strip().split("\t")[2:]
	for s in samplesList[2:]:
		newCount=0
		#print s
		samp=s.strip()+"_"+str(samplesList.index(s))
		origCount=float(counts[samplesIndexLookup[samp]])
		#print origCount
		newCount=float(origCount)*scalingFactors[samp]
		# check output that these values match column sums (col sums are correct)
		if samp not in newLibSizes:
			#print samp, "x"
			newLibSizes[samp]=0
			newLibSizes[samp]+=float(newCount)
			# print newLibSizes
		else:
			newLibSizes[samp]+=float(newCount)
		OFH.write(str(newCount))
		OFH.write("\t")
	OFH.write("\n")

print "new lib sizes"
for s in sorted(newLibSizes):
	print s, newLibSizes[s]

f.close()
OFH.close()
			
		
