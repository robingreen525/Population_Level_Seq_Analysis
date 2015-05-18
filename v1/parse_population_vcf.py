import sys

f=open(sys.argv[1],'r')

ll=f.readlines()

for line in ll:
	line=line.split()
        #print(line)
	if line[0][0]=='S' or line[0][0]=='c':
		#print line[0],line[1],line[3],line[4]
		info=line[7].split(';')
                #print(info)
		#print info[0],info[5],info[7]
		prop=info[0].split('=')[1]
		print line[0],line[1],line[3],line[4],prop
