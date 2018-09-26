#!/usr/bin/env python

def pullseqs(GIfile):
	from Bio import Entrez, SeqIO
	import os, sys, time

	Entrez.email = "sunray1@ufl.edu"
	seqids = []
	records = []

	name = GIfile.split(".")[0]
	with open(GIfile) as o:
		line = o.readline()
		while line:
			seqids.append(line.strip())
			line = o.readline()
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
		print(str(round((float(len(records))/float(len(seqids)))*100, 2)) + "%")
		handle.close()

	SeqIO.write(records, name + ".fa", "fasta")

