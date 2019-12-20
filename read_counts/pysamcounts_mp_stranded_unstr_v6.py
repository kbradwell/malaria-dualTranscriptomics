import sys
from sys import argv
import re
import pysam
from collections import Counter
import glob
import multiprocessing

# outputs counts for each gene of all reads covering any of its exons, in any of the isoforms
# v6: changes the loops for unstr and stranded input so unnecessary calcs aren't performed

# MULTIPROCESSING
# split a list into evenly sized chunks
def chunks(l, n):
    return [l[i:i+n] for i in range(0, len(l), n)]

def do_job(job_id, data_slice, list_of_files, gftHash, chrList, samplesList, libType):
	#print "job", job_id, data_slice,
	#print "extracting unique reads"
	geneReadsHash={}
	seen={}
	for gene in data_slice:
		for file_name in list_of_files:
			samfile = pysam.Samfile(file_name, "rb")
			if str(file_name) not in geneReadsHash:
				geneReadsHash[str(file_name)]={}
			if str(file_name) not in seen:
				seen[str(file_name)]={}
			for exonTup in gftHash[gene]:
				for read in samfile.fetch(exonTup[0],int(exonTup[2])-1, int(exonTup[3])+1):
					if read.qname not in seen[str(file_name)]:
						#print dir(read)
						if "unstr" in libType:
							if gene not in geneReadsHash[str(file_name)]:
								geneReadsHash[str(file_name)][gene]=1
							else:
								geneReadsHash[str(file_name)][gene]+=1
						else: # it's stranded data
							NumInPair=0
							binFlag = ('{0:b}'.format(int(read.flag)))[::-1]
							flagArray=list(binFlag)
							if flagArray[6]=="1":
								NumInPair=1
								if "+" in gftHash[gene][0][1]:
									if flagArray[4]=="0":
										if gene not in geneReadsHash[str(file_name)]:
											geneReadsHash[str(file_name)][gene]=1
										else:
											geneReadsHash[str(file_name)][gene]+=1
								elif "-" in gftHash[gene][0][1]:
									if flagArray[4]=="1":
										if gene not in geneReadsHash[str(file_name)]:
											geneReadsHash[str(file_name)][gene]=1
										else:
											geneReadsHash[str(file_name)][gene]+=1
							elif flagArray[7]=="1":
								NumInPair=2
								if "+" in gftHash[gene][0][1]:
									if flagArray[4]=="1":
										if gene not in geneReadsHash[str(file_name)]:
											geneReadsHash[str(file_name)][gene]=1
										else:
											geneReadsHash[str(file_name)][gene]+=1
								elif "-" in gftHash[gene][0][1]:
									if flagArray[4]=="0":
										if gene not in geneReadsHash[str(file_name)]:
											geneReadsHash[str(file_name)][gene]=1
										else:
											geneReadsHash[str(file_name)][gene]+=1
						seen[str(file_name)][read.qname]=1
							

	#print "printing counts"
	for g in data_slice:
		print "\n",g,"\t","-",
		for s in sorted(samplesList):
			print "\t",
			if s in geneReadsHash:
				if g in geneReadsHash[s]:
					print str(geneReadsHash[s][g]),
				else:
					print "0",
				
def dispatch_jobs(data, job_number, infiles, coords, chrList, samples, libType):
	total = len(data)
	chunk_size = total / job_number
	slice = chunks(data, chunk_size)
	jobs = []

	for i, s in enumerate(slice):
		j = multiprocessing.Process(target=do_job, args=(i, s, infiles, coords, chrList, samples, libType))
		jobs.append(j)
	for j in jobs:
		j.start()
	for j in jobs:
		j.join()


if __name__ == "__main__":

	script, gtffile, libType = argv

	GTF = open(gtffile,'r')

	# read in exon coords for each gene from GTF file
	testflag =('{0:b}'.format(131))[::-1]
	testflagArray=list(testflag)
	#print testflagArray
	gftHash={}
	geneNameList=[]
	chrList=[]
	for line in GTF:
		coordParts=line.strip().split("\t")
		chr=coordParts[0]
		if chr not in chrList:
			chrList.append(chr)
		featureType=coordParts[2]
		strand=coordParts[6]
		if "exon" in featureType:
			start=coordParts[3]
			end=coordParts[4]
			geneNameAnnot=coordParts[8].strip()
			geneName=geneNameAnnot.split(";")[1].split("\"")[1]
			#print geneName
			if geneName not in geneNameList:
				geneNameList.append(geneName)
			if geneName not in gftHash:
				gftHash[geneName]=[]
				coordTup=(chr,strand,int(start),int(end))
				gftHash[geneName].append(coordTup)
			else:
				coordTup=(chr,strand,int(start),int(end))
				gftHash[geneName].append(coordTup)

	list_of_files = glob.glob('./*nodup.bam')  # create the list of files

	samplesList=[]
	for file_name in list_of_files:
		samplesList.append(str(file_name))
			
	print "GeneID\tLength",
	for s in sorted(samplesList):
			print "\t",s,

	data = geneNameList
	infiles = list_of_files
	coords = gftHash
	samples = samplesList
	if len(data)>500:
		geneChunks=chunks(data,100)
		#print geneChunks
		for g in geneChunks:
			coordChunk = { g_key: coords[g_key] for g_key in g }
			#print coordChunk.values()
			chrChunk=[]
			for x in coordChunk.values():
				geneChr=[locTup[0] for locTup in x]
				for chr in geneChr:
					chrChunk.append(chr)
			chrChunkSet=set(chrChunk)
			#print "dispatching jobs"
			dispatch_jobs(g, 5, infiles, coordChunk, chrChunkSet, samples, libType)
	else:
		dispatch_jobs(data, 5, infiles, coords, chrList, samples, libType)
	GTF.close()
