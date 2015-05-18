import sys

f=open(sys.argv[1],'r')

rm11_genes={}

ll=f.readlines()
for line in ll:
	line=line.split('\t')
	#print line[0],line[1]
	if line[0] not in rm11_genes.keys():
			if line[1]=='':
				rm11_genes[line[0]]='Unknown'
			else:
				
				rm11_genes[line[0]]=str(line[1])




f=open(sys.argv[2],'r')
ll=f.readlines()

contig_sites={}

for line in ll:
	#line=ll[i]
        #line=line.rstrip('\n')
	line=line.split('\t')
        print line
	gene_id=line[8].split(';')[0].split()[1]
	if line[2]=='exon':
		start=int(line[3])
		end=int(line[4])
		contig=line[0]
		#print gene_id
		if contig not in contig_sites.keys():
				contig_sites[contig]=[start,end,gene_id]
		else:
				contig_sites[contig].append(start)  # making a dictionary of list [star1,end1, star2,end2,....] for each contig
				contig_sites[contig].append(end)
				contig_sites[contig].append(gene_id)


f=open(sys.argv[3],'r')
m=open(sys.argv[4],'r')

auxo={}

ll=m.readlines()
for line in ll:
	line=line.split('\t')
	if len(line)>1:
		auxo[line[0]]=line[7]









ll=f.readlines()
for line in ll:
		line1=line
		line=line.split()
		contig=line[0]
		site=int(line[1])
		#print contig,site
		i=0
		while i <(len(contig_sites[contig])-2):
				start=contig_sites[contig][i]
				end=contig_sites[contig][i+1]
				id=contig_sites[contig][i+2]
				
				if site>=start and site <=end:
						id= id.split('"')[1]
						#print(len(rm11_genes[id]))
						#print rm11_genes[id]
						temp=[]
						for item in line:
							if item !='':
								temp.append(item.rstrip('\r\n'))
						
						#print rm11_genes[id]
						temp.append(rm11_genes[id])
						gene=rm11_genes[id]
						if gene in auxo.keys():
							temp.append(auxo[gene])

						else:
							temp.append('Not Auxotroph')
						
						print temp

				i=i+3



f.close()
