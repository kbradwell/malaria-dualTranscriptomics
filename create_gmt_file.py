# create custom gmt file

PLASMODIUM_GENES_IFH = open("./plasmodium_genes_list.txt",'r')
UNIPROT_PLASMODIUMID_IFH = open("./uniprot_plasmodbID.txt",'r')
UNIPRO_PATHWAYID_IFH = open("./plasmodium_uniprot_id_to_reactome_pathway.txt",'r')
PATHWAYID_PATHWAYNAME_IFH = open("./plasmodium_list_pathways_names_reactome.txt",'r')
OFH = open("plasmodium.gmt",'w')

# format of gmt file:
# pathway name \t description (website) \t tab sep gene list

plasmodium_genes = []

for l in PLASMODIUM_GENES_IFH:
	gene = l.strip()
	plasmodium_genes.append(gene)

mappedPlasmoID = []
uniprotHash = {}
for l in UNIPROT_PLASMODIUMID_IFH:
	genemap = l.strip()
	cols = genemap.split('\t')
	uniprotID = cols[0]
	plasmodiumID = cols[1]
	mappedPlasmoID.append(plasmodiumID)
	if uniprotID not in uniprotHash:
		uniprotHash[uniprotID] = [plasmodiumID]
	else:
		uniprotHash[uniprotID].append(plasmodiumID)

print(len(set(mappedPlasmoID)))
print(len(set(plasmodium_genes)))

pathwayIDhash = {}
for l in UNIPRO_PATHWAYID_IFH:
	idmap = l.strip()
	cols = idmap.split('\t')
	uniprotID = cols[0]
	pathwayID = cols[1]
	if pathwayID not in pathwayIDhash:
		pathwayIDhash[pathwayID] = [uniprotID]
	else:
		pathwayIDhash[pathwayID].append(uniprotID)

pathwayID2Name = {}
for l in PATHWAYID_PATHWAYNAME_IFH:
	namemap = l.strip()
	cols = namemap.split('\t')
	pathwayID = cols[0]
	pathwayName = cols[1]
	if pathwayID not in pathwayID2Name:
		pathwayID2Name[pathwayID] = pathwayName

for pid in pathwayID2Name:
	if pid in pathwayIDhash:
		OFH.write(pathwayID2Name[pid])
		OFH.write("\t")
		OFH.write(".")
		uniprotIDlist = pathwayIDhash[pid]
		for uniprotID in uniprotIDlist:
			for plasmodbID in uniprotHash[uniprotID]:
				OFH.write("\t")
				OFH.write(plasmodbID)
		OFH.write("\n")
	
	