import sys

f=open(sys.argv[1],'r')

ll=f.readlines()


genes=[]

for line in ll:
	line=line.split('\t')
        #print line[0]
	if line[0]=='NC_012967.1':
		#print line
                if line[2]=='gene':
                        #print line
                        temp=[line[0],int(line[3]),int(line[4])]
                        #print temp
                        alias=line[8].split(';')
                        #print alias
                        for entry in alias:
                                entry=entry.split('=')
                                if entry[0]=='Name':
                                        gene=entry[1]
                                        temp.append(gene)
                                        genes.append(temp)

f.close()

print('herp')

f=open(sys.argv[2],'r')

ll=f.readlines()

#print ll

for line in ll:
        line=line.split('\t')
        #print line
        if line[0][0:2]=='gi':
                chr=line[0]
                #print line
                #print chr
                site=int(line[1])
                print site
                for gene in genes:
                        if site>=gene[1]:
                                if site <= gene[2]:
                                        #print gene
                                        print gene[3]
                                        print line,gene[3]
