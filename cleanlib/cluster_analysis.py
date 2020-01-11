#!/usr/bin/env python
#does cluster analysis for those I couldn't resolve by dna content or length
import sys, sqlite3, time
from Bio import SeqIO, Entrez, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.Align.Applications import MuscleCommandline
from Bio.Align import AlignInfo
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
from blastlib.clean_seq_funcs import alignment_comp, alignment_reg, alignment_rev_comp, blast, identity_calc, tiling
from cleanlib.databasing import get_seqs_from_sqldb_GI_no_gene, export_fasta, local_blast

def cluster(blastdb, taxdb, email):
    Entrez.email = email
    conn = sqlite3.connect(blastdb)
    c = conn.cursor()
    c.execute("ATTACH '" + taxdb + "' as 'tax'")
    muscle_cline = MuscleCommandline(clwstrict=True)
    input_dic = {}
    multiple_dic = {}
    two_dic = {}
    problem_dic = {}
    finalseqs = set()
    multfinalseqs = []
    unresolved = []
    
    for iter in c.execute("SELECT tc_id, GI from blast where decision = 'To cluster analysis/Chosen' order by tc_id"):
        if iter[0] not in input_dic:
            input_dic[iter[0]]  = [str(iter[1])]
        else:
            output = input_dic[iter[0]]
            output.append(str(iter[1]))
            input_dic[iter[0]] = output
        
    for i in input_dic:
        GIs_list = input_dic[i]
        if len(GIs_list) > 2:
            multiple_dic[i] = GIs_list
        if len(GIs_list) == 2:
            two_dic[i] = GIs_list
    
    
    for i in multiple_dic:
        identities = []
        #get consensus of all
        for iter in c.execute("SELECT Gene_name FROM blast WHERE GI = '"+ str(multiple_dic[i][0]) +"';"):
            gene = str(iter[0])
        align = alignment_reg(multiple_dic[i], blastdb, False, gene, c)
        summary_align = AlignInfo.SummaryInfo(align)
        consensus = summary_align.gap_consensus(threshold=.5, ambiguous='N')
        consensus_record = SeqRecord(consensus, id="Consensus_all")
        for m in multiple_dic[i]:
            #align each gene to consensus
            seqs = []
            iterator = get_seqs_from_sqldb_GI_no_gene([m], "hseq", blastdb, c)
            for seq in iterator:
                seqs.append(seq)
            handle_string = StringIO()
            SeqIO.write(seqs, handle_string, "fasta")
            SeqIO.write(consensus_record, handle_string, "fasta")
            data = handle_string.getvalue()
            stdout, stderr = muscle_cline(stdin=data)
            align = AlignIO.read(StringIO(stdout), "clustal")
            count = 0
            gaps = 0
            for col in range(0, len(align[0])):
                column = align[:,col]
                if "-" not in column:
                    if column[1:]==column[:-1]:
                        count=count+1
                else:
                    gaps=gaps+1
            iden = 100*(count/float((len(align[0])-gaps)))
            identities.append(iden)
        if identities.count(max(identities)) == 1:
            GI_to_pick = multiple_dic[i][identities.index(max(identities))]
            c.execute("UPDATE blast SET Decision='Closest to consensus in cluster analysis/Chosen' WHERE GI = '" + GI_to_pick + "';")
            GIS_not_picked = list(set(multiple_dic[i])-set([GI_to_pick]))
            GIS_not_picked_str = str(GIS_not_picked).replace("[", "(").replace("]", ")")
            c.execute("UPDATE blast SET Decision='Further from consensus in cluster analysis/Not chosen' WHERE GI IN " + GIS_not_picked_str + ";") 
        else:
            GI_to_pick = [multiple_dic[i][m] for m, x in enumerate(identities) if x == max(identities)]
            GI_to_pick_str = str(GI_to_pick).replace("[", "(").replace("]", ")")
            c.execute("UPDATE blast SET Decision='Pick one randomly/Chosen' WHERE GI IN " + GI_to_pick_str + ";")
            
    for i in two_dic:
        twodic_str = str(two_dic[i]).replace("[", "(").replace("]", ")")
        c.execute("UPDATE blast SET Decision='Pick one randomly/Chosen' WHERE GI IN " + twodic_str + ";")

    conn.commit()
    conn.close()
