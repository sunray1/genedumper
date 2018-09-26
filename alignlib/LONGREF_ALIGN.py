#! /usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Alphabet.IUPAC import ambiguous_dna, extended_protein

import numpy
import sys
import os


def align(fastain):
#	thread = arguments[arguments.index("-T")+1]
	thread = 2
	outparts=fastain.split(".")
	seqL=[]
	runfull = "mafft --adjustdirection --maxiterate 1000 --localpair --thread %s %s > %s"
	addfrag = "mafft --adjustdirection --maxiterate 1000 --localpair --thread %s --addfragments %s %s > %s"
	killfile = "rm %s"
	
	fullcount=0
	fragcount=0
	
	
	fasta_sequences = SeqIO.parse(open(fastain),'fasta')
	end = False
	outfile1 =outparts[0] + "_FULL.fas"
	outfile2 =outparts[0] + "_FRAGS.fas"
	outfile3 =outparts[0] + "_FULLALIGN.fa"
	outfile4 =outparts[0] + "_FINALALIGN.fa"
	
	for record in fasta_sequences:
		record.description=""
		SEQ = record.upper()
		seqL.append(len(SEQ.seq))
	
	#print seqL
	dev=numpy.std(seqL)*2
	upper=max(seqL)
	above=upper-dev
	
	f= open(outfile1, "w")
	g= open(outfile2, "w")
	
	fasta_sequences = SeqIO.parse(open(fastain),'fasta')
	end = False
	for record in fasta_sequences:
		record.description=""
		SEQ = record.upper()	
		if len(SEQ.seq) >= above:
			SeqIO.write([record], f, "fasta")		
			fullcount+=1
		else:
			SeqIO.write([record], g, "fasta")
			fragcount+=1
	
	f.close()
	g.close()
	
	print(outfile1 + ":" + str(fullcount) + " seqs\n" + outfile2 + ": " + str(fragcount) + " seqs\n")
	print("RUNNING MAFFT FULL ALIGNMENT.....")
	os.system(runfull % (thread, outfile1, outfile3))
	print("RUNNING MAFFT ADDING FRAGS.....")
	os.system(addfrag % (thread, outfile2, outfile3, outfile4))
	os.system(killfile %(outfile1))
	os.system(killfile %(outfile2))
	os.system(killfile %(outfile3))
