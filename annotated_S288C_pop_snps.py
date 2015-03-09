import sys

f=open(sys.argv[1],'r')

ll=f.readlines()


genes=[]

for line in ll:
	line=line.split('\t')
	if line[0][0:3]=='chr':
		#print line
                if line[2]=='gene':
                        #print line
                        temp=[line[0],int(line[3]),int(line[4])]
                        #print temp
                        alias=line[8].split(';')
                        for entry in alias:
                                entry=entry.split('=')
                                if entry[0]=='gene':
                                        gene=entry[1]
                                        temp.append(gene)
                                        genes.append(temp)

f.close()

f=open(sys.argv[2],'r')

ll=f.readlines()

for line in ll:
        line=line.split('\t')
        if line[0][0:3]=='chr':
                chr=line[0]
                site=int(line[1])
                for gene in genes:
                        if chr==gene[0]:
                                if site>=gene[1]:
                                        if site <= gene[2]:
                                                #print gene
                                                print gene[3]
                                                print line,gene[3]
