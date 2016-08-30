###############
# 2015-07-05
# the purpose of this script is to download (in bulk) the Refseq IDs 
# for all yeast proteins in the NCBI Entrez database.
# This script uses built in functionalities from the BioPython libraries
# I will use this information to submit my non-synonymous mutations to SIFT Protein in 
# bulk to try and get predictions on whether these mutations may be loss of function.

# usuage: python get_refseq_yeast.py

# Note: using the sbatch grab_RefSeq.sh command to run this script on the cluster

# Biopython Manual: http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc130
# Entrez Package Documentation: http://biopython.org/DIST/docs/api/Bio.Entrez-module.html

# import Biopython libraies
from Bio import Entrez

# information for my reference
entrezDbName = 'protein'
ncbiTaxId = '4932' # S. cerevisea 
ncbiTaxId= '559292' # S. cerevisiea s288C
Entrez.email = 'rjgreen@fhcrc.org'


# Find entries matching the query
entrezQuery = "refseq[filter] AND txid559292[Organism]" # actual query 
searchResultHandle = Entrez.esearch(db=entrezDbName, term=entrezQuery,usehistory='y') # run query against Entrez db, use history
# allows for running things in small processes
searchResult = Entrez.read(searchResultHandle) # Parses an XML file from the NCBI Entrez Utilities into python objects.
searchResultHandle.close()


# these parameters are needed for getting sequences in small batches of 100
webenv = searchResult["WebEnv"]
query_key= searchResult['QueryKey']
count = int(searchResult["Count"])

#print(searchResult)

# loop for get all proteins seq, 100 at a time (see section 9.15 of BioPython Manual, this is the exact same setup as example) 
batch_size = 10
out_handle = open("S.cerevisea_RefSeq.fasta", "w") # save sequences to this fasta file
for start in range(0,count,batch_size):
	end = min(count, start+batch_size)
	print("Going to download record %i to %i" % (start+1, end))
	fetch_handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text",
                                 retstart=start, retmax=batch_size,
                                 webenv=webenv, query_key=query_key)
	data = fetch_handle.read()
	fetch_handle.close()
	out_handle.write(data)
    
out_handle.close()

