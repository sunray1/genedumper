#!/usr/bin/env python

def pullseqs(blastdb, email):
	from Bio import Entrez, SeqIO
	import os, sys, time, sqlite3
	from cleanlib.databasing import get_seqs_from_sqldb_GI
	conn = sqlite3.connect(blastdb)
	c = conn.cursor()
	GI_gene_dic = {}
	mitochloro_gene_dic = {}
	genes = set()
	Entrez.email = email
# Ambiguous species/not chosen
# Better tiling/Not chosen
# Chosen, but does not reciprocal blast
# Chosen, but does not reciprocal blast and all do not align
# Closest to consensus in cluster analysis/Chosen
# Further from consensus in cluster analysis/Not chosen
# Longest or most info, good top hit/chosen
# Mito or chloro sequence/Chosen
# Only choice/chosen
# Pick one randomly/Chosen
# Sequence did not have same top blast species, but all aligned correctly/Chosen
# Short or less info/Not chosen
# Species not in taxonomy/not chosen
# To cluster analysis/Chosen
	
	for iter in c.execute("SELECT GI, Gene_name FROM blast WHERE Decision IN ('Closest to consensus in cluster analysis/Chosen', 'Longest or most info, good top hit/chosen', 'Only choice/chosen', 'Sequence did not have same top blast species, but all aligned correctly/Chosen')"):
		GI_gene_dic[str(iter[0])] = iter[1]
		genes.add(iter[1])
	
	for iter in c.execute("SELECT GI, Gene_name FROM blast WHERE Decision IN ('Mito or chloro sequence/Chosen')"):
		mitochloro_gene_dic[str(iter[0])] = iter[1]	
		genes.add(iter[1])
	
		
	for gene in genes:
		#get regular sequences
		records = []
		seqids = [i for i in GI_gene_dic if GI_gene_dic[i] == gene]
		seqlists = [seqids[i:i+200] for i in range(0, len(seqids), 200)]
		for i in seqlists:
			error = False
			seqids_sub = ",".join(i)
			try:
				handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=seqids_sub)
			except:
				error = True
			while error is True:
				try:
					print("Error, trying again")
					time.sleep(5)
					handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=seqids_sub)
					error = False
				except:
					pass
			for seq_record in SeqIO.parse(handle, "fasta"):
				records.append(seq_record)
				
		GI_mito_GI = [i for i in mitochloro_gene_dic if mitochloro_gene_dic[i] == gene]
		iterator = get_seqs_from_sqldb_GI(GI_mito_GI, "hseq", blastdb, gene, c)
		for seq in iterator:
			records.append(seq)
		# print(str(round((float(len(records))/float(len(seqids)))*100, 2)) + "%")


		SeqIO.write(records, gene + "_cleaned.fa", "fasta")
	conn.close()


