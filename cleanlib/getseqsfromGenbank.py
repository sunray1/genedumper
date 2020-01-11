#!/usr/bin/env python

def pullseqs(blastdb, email):
	from Bio import Entrez, SeqIO
	import os, sys, time, sqlite3
	from random import randint
	from cleanlib.databasing import get_seqs_from_sqldb_GI
	conn = sqlite3.connect(blastdb)
	c = conn.cursor()
	GI_gene_dic = {}
	mitochloro_gene_dic = {}
	tc_ids_random = set()
	genes = set()
	Entrez.email = email

# Ambiguous species/not chosen
# Better tiling/Not chosen
# Closest to consensus in cluster analysis/Chosen
# Further from consensus in cluster analysis/Not chosen
# Longest or most info, good top hit/chosen
# Mito or chloro sequence/Chosen
# Only choice/chosen
# Only or best choice in tiling analysis/chosen
# Pick one randomly/Chosen
# Sequence did not have same top blast species, but all aligned correctly, tile 1/Chosen
# Short or less info, tile 1/Not chosen
# Short or less info/Not chosen
# Species not in taxonomy/not chosen
	
	for iter in c.execute("SELECT GI, Gene_name FROM blast WHERE Decision IN ('Closest to consensus in cluster analysis/Chosen', 'Longest or most info, good top hit/chosen', 'Only choice/chosen', 'Only or best choice in tiling analysis/chosen') OR Decision LIKE 'Sequence did not have same top blast species, but all aligned correctly%'"):
		GI_gene_dic[str(iter[0])] = iter[1]
		genes.add(iter[1])
	
	for iter in c.execute("SELECT GI, Gene_name FROM blast WHERE Decision IN ('Mito or chloro sequence/Chosen')"):
		mitochloro_gene_dic[str(iter[0])] = iter[1]	
		genes.add(iter[1])
	
		
	for gene in genes:
		#get regular sequences
		pick_random_dic = {}
		records = []
		#get regular
		
		seqids = [i for i in GI_gene_dic if GI_gene_dic[i] == gene]
		
		#choose random
		for iter in c.execute("SELECT tc_id, GI FROM blast WHERE Decision IN ('Pick one randomly/Chosen') AND Gene_name = '"+gene+"' ORDER BY tc_id"):
			if iter[0] not in pick_random_dic:
				pick_random_dic[iter[0]] = [iter[1]]
			else:
				output = pick_random_dic[iter[0]]
				output.append(iter[1])
				pick_random_dic[iter[0]] = output
		
		for i in pick_random_dic:
			choice = randint(0, len(pick_random_dic[i])-1)
			GI_choice = pick_random_dic[i][choice]
			seqids.append(str(GI_choice))
			
			
		#pull the seqs
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
		#get mito/chloro		
		GI_mito_GI = [i for i in mitochloro_gene_dic if mitochloro_gene_dic[i] == gene]
		if len(GI_mito_GI) != 0:
			iterator = get_seqs_from_sqldb_GI(GI_mito_GI, "hseq", blastdb, gene, c)
			for seq in iterator:
				records.append(seq)
		
			
		# print(str(round((float(len(records))/float(len(seqids)))*100, 2)) + "%")


		SeqIO.write(records, gene + "_cleaned.fa", "fasta")
	conn.close()


