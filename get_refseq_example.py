from Bio import Entrez

entrezDbName = 'protein'
#ncbiTaxId = '1001533' # Bovine papillomavirus 7
ncbiTaxId = '4932' # S. cerevisea 

Entrez.email = 'rjgreen@fhcrc.org'

# Find entries matching the query
entrezQuery = "refseq[filter] AND txid%s"%(ncbiTaxId)
searchResultHandle = Entrez.esearch(db=entrezDbName, term=entrezQuery)
searchResult = Entrez.read(searchResultHandle)
searchResultHandle.close()

# Get the data.
uidList = ','.join(searchResult['IdList'])
entryData = Entrez.efetch(db=entrezDbName, id=uidList, rettype='fasta').read()
print entryData