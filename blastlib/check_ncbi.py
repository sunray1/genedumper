#!/usr/bin/env python
#Script to run step 4 - designed to check ncbi to make sure the names parsed out of the blast hits are the names listed as
# the organism names in the nucleotide database. The species name in the blast database will be updated if the ncbi organism name
# is in the taxonomy database.
#
# Changes are noted in 'ncbi.txt'
def ncbi(taxdb, blastdb, email):
	import time, sqlite3, sys
	from Bio import Entrez
	Entrez.email = email
	num_spe_dic = {}
	count = 1
	conn = sqlite3.connect(taxdb)
	c = conn.cursor()
	c.execute("ATTACH '" + blastdb + "' as 'db'")
	#make dictionary of species and accession numbers for all records where the nc_id is null
	for iter in c.execute("SELECT accession, Species FROM blast WHERE tc_id IS NULL;"):
		num_spe_dic[str(iter[0])] = str(iter[1])
	with open('ncbi.txt', 'w') as o:
		chunksncbi = [num_spe_dic.keys()[x:x+5000] for x in range(0, len(num_spe_dic.keys()), 5000)]
		print(len(chunksncbi))
		for chunk in chunksncbi:
			print(count)
			ncbin = ",".join(chunk)
			try:
				handle = Entrez.efetch(db='nucleotide', id=ncbin, retmode='xml', seq_start=1, seq_stop=1)
				print("parsing records")
				records = Entrez.read(handle)
				print("finished parsing")
			except:
				print("killed early, trying again")
				handle = Entrez.efetch(db='nucleotide', id=ncbin, retmode='xml', seq_start=1, seq_stop=1)
				print("parsing records")
				records = Entrez.read(handle)
				print("finished parsing")				
			for r in records:
				namefind = 'None'
				num = r[u'GBSeq_accession-version']
				organism = r[u'GBSeq_organism'].replace("'", "")
				if num_spe_dic[num] == organism:
					pass
				else:
		#if the genbank name and original name are not the same AND the genbank name is in our list - change original name
			#to genbank name in sqlite database
					for iter in c.execute("SELECT namestr FROM names WHERE namestr = '" + organism + "';"):
						namefind = iter[0]
				if namefind != 'None':
					print('updating records')
					c.execute("UPDATE blast SET Species='" + organism + "' WHERE accession='" + num + "';")
					c.execute("UPDATE blast SET genus='" + organism.split()[0] + "' WHERE accession='" + num + "';")
					c.execute("UPDATE blast SET epithet='" + organism.split()[1] + "' WHERE accession='" + num + "';")
					o.write(num_spe_dic[num] + '\t' + num + '\t' + organism + '\n')
				handle.close()
			conn.commit()
			count += 1
	conn.commit()
	conn.close()		
