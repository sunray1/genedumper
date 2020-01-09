#!/usr/bin/env python
#script designed to create and execute local blast for cleaning

def get_seqs_from_sqldb_GI(GIset, seqtype, blastdb, gene, c):
	#returns iterator of seqrecs
    #seqtype = qseq or hseq
    import os, sys, sqlite3
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    GIset = str(GIset).replace("[", "(").replace("]", ")")
    if seqtype == "hseq":
        sql = "SELECT accession, "+ seqtype +" FROM blast WHERE Gene_name = '" + gene + "' AND GI IN "+ GIset +";"
        for iter in c.execute(sql):
            record = SeqRecord(Seq(str(iter[1])), id = iter[0],  description = "")
            yield(record)
    if seqtype == "qseq":
        for iter in c.execute("SELECT Gene_name, "+ seqtype +" FROM blast WHERE GI IN "+ GIset +";"):
            record = SeqRecord(Seq(str(iter[1])), id = iter[0],  description = "")
            yield(record)
            
def get_seqs_from_sqldb_GI_no_gene(GIset, seqtype, blastdb, c):
	#returns iterator of seqrecs
    #seqtype = qseq or hseq
    import os, sys, sqlite3
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    GIset = str(GIset).replace("[", "(").replace("]", ")")
    if seqtype == "hseq":
        sql = "SELECT accession, "+ seqtype +" FROM blast WHERE GI IN "+ GIset +";"
        for iter in c.execute(sql):
            record = SeqRecord(Seq(str(iter[1])), id = iter[0],  description = "")
            yield(record)
    if seqtype == "qseq":
        for iter in c.execute("SELECT Gene_name, "+ seqtype +" FROM blast WHERE GI IN "+ GIset +";"):
            record = SeqRecord(Seq(str(iter[1])), id = iter[0],  description = "")
            yield(record)



def get_seqs_from_sqldb(accset, seqtype, blastdb, c):
	#returns iterator of seqrecs
    #seqtype = qseq or hseq
    import os, sys, sqlite3
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    accset = str(list(accset)).replace("{", "(").replace("}", ")").replace("[", "(").replace("]", ")")
    for iter in c.execute("SELECT accession, "+ seqtype +" FROM blast WHERE accession IN "+ accset +";"):
        record = SeqRecord(Seq(str(iter[1])), id = iter[0],  description = "")
        yield(record)
#usage
#seqs = get_seqs_from_sqldb(["1080121667", "914615309"], "hseq", "blast_results.db")
#print(next(seqs))
#print(next(seqs))

def export_fasta(Seqiterator, outfilename):
    from Bio import SeqIO
    SeqIO.write(Seqiterator, outfilename, "fasta")
    
    
def create_blast_db(fastafile):
    from Bio.Blast.Applications import NcbimakeblastdbCommandline
    cline = NcbimakeblastdbCommandline(dbtype="nucl", input_file=fastafile)
    stdout, stderr = cline()

def local_blast(dbfile, queryfile):
    from Bio.Blast.Applications import NcbiblastnCommandline
    outfile = queryfile+".xml"
    blastn_cline = NcbiblastnCommandline(query=queryfile, db=dbfile, evalue = 0.001, word_size=28, outfmt=5, out=outfile)
    stdout, stderr = blastn_cline()


