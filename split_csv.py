from sys import argv
import os

script, csvfile = argv

IFH=open(csvfile,'r')

linesHg=12487

header = IFH.readline()

linecounter = 0

hgLines=[]
pfLines=[]

for line in IFH:
	if linecounter<linesHg:
		hgLines.append(line)
	else:
		pfLines.append(line)
	linecounter+=1
	
	
OFH1=open('split_pf.csv','w')

for elem in pfLines:
	OFH1.write(elem)
	
OFH1.close()

n=781 
LoL=[hgLines[i:i + n] for i in xrange(0, len(hgLines), n)]

filecount = 1
for chunk in LoL:
	outfilebase = os.path.basename(csvfile)
	outfilename = str(outfilebase).split(".csv")[0]+"_hggenes_file"+str(filecount)+"_split.csv"
	OFH=open(outfilename,'w')
	for elem in chunk:
		OFH.write(elem)
	filecount+=1
	OFH.close()

IFH.close()
