import sys

cutoff=float(sys.argv[3])
f=open(sys.argv[1],'r')

ll=f.readlines()

auxo={}


##build database of known auxotrophy mutations from SGD 2015 dataset (auxotrophy_annotations.txt)
for line in ll:
    line=line.split('\t')
    if len(line)>5:
        #print line[0],line[7]
        if line[0] not in auxo.keys():
            auxo[line[0]]=line[7]


f.close()

#print auxo.keys()

# find auxotroph snps >= set cutoff
f=open(sys.argv[2],'r')

ll=f.readlines()

for line in ll:
    line=line.rstrip('\n')
    line1=line
    line =line.split(',')
    gene=line[len(line)-1].split(']')
    if len(gene)>1:
        gene=gene[1][1:]
        if gene in auxo.keys():
           # print gene,auxo[gene]
            line=line1.split(',')[7]
            if len(line)>100: ## make sure this line isn't causing me to lose information
                freq=line.split(';')[0]
                alt_reads=line.split(';')[5]
                all_reads=line.split(';')[28]
                freq=freq.split('=')
                if len(freq)>1:
                    freq=float(freq[1])
                    if freq>=cutoff:
                        print gene,auxo[gene],freq,alt_reads,all_reads

                    
