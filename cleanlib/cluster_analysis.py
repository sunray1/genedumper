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
from cleanlib.databasing import get_seqs_from_sqldb_GI, export_fasta, local_blast

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
    
    with open("multiple_gene_choices.txt") as o:
        line = o.readline()
        while line:
            input_dic[line.split("\t")[0]] = line.strip().split("\t")[1].replace("[", "").replace("]", "").replace("'", "")
            line = o.readline()
    for i in input_dic:
        GIs = input_dic[i]
        GIs_list = GIs.split(", ")
        if len(GIs_list) > 2:
            multiple_dic[i] = GIs_list
        if len(GIs_list) == 2:
            two_dic[i] = GIs_list
    
    
    for i in multiple_dic:
        identities = []
        #get consensus of all
        for iter in c.execute("SELECT Gene_name FROM blast WHERE GI = '"+ multiple_dic[i][0] +"';"):
            gene = str(iter[0])
        align = alignment_reg(multiple[i], blastdb, "False", gene, c)
        summary_align = AlignInfo.SummaryInfo(align)
        consensus = summary_align.gap_consensus(threshold=.5, ambiguous='N')
        consensus_record = SeqRecord(consensus, id="Consensus_all")
        for m in multiple_dic[i]:
            #align each gene to consensus
            iterator = get_seqs_from_sqldb_GI(m, "hseq", blastdb, gene, c)
            handle_string = StringIO()
            for seq in iterator:
                seqs.append(seq)
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
            finalseqs.add(multiple_dic[i][identities.index(max(identities))])
        else:
            problem_dic[i] = multiple_dic[i]
            GI_to_pick = [multiple_dic[i][m] for m, x in enumerate(identities) if x == max(identities)]
            multfinalseqs.append(GI_to_pick)
        
    for i in two_dic:
        #align the two seqs
        list_of_GIs = two_dic[i]
        for iter in c.execute("SELECT Gene_name FROM blast WHERE GI = '"+ list_of_GIs[0] +"';"):
            gene = str(iter[0])
        alignment = alignment_reg(list_of_GIs, blastdb, "False", gene, c)
        iden = identity_calc(alignment)
        if iden < 95:
#            print("Low Aligned Identity: " + str(iden))
            alignment = alignment_rev_comp(list_of_GIs, blastdb, "False", gene, c)
            iden = identity_calc(alignment)
            if iden < 95: 
    #get taxonomy for query(main species)
 #               print("Low Reverse Complement Aligned Identity: " + str(iden))
                alignment = alignment_comp(list_of_GIs, blastdb, "False", gene, c)
                iden = identity_calc(alignment)
                if iden < 95:
 #                   print("Low Complement Aligned Identity: " + str(iden))
            #add tiling thing
                    gene_name = '_'.join(i.split('_')[1:])
                    idens, start_stop = tiling(list_of_GIs, gene_name, blastdb, c)
                    current_start = -1
                    current_stop = -1 
                    result = []
                    if all(m > 70 for m in idens):
                        for start, stop in sorted(start_stop):
                            if start > current_stop:
                                result.append((start, stop))
                                current_start, current_stop = start, stop
                            else:
                                current_stop = max(current_stop, stop)
                                result[-1] = (current_start, current_stop)
                        if len(result) == len(start_stop):
 #                           print("Seqs align to different regions of probe, choosing all")
                            multfinalseqs.append(list_of_GIs)
                        else:
 #                           print('Seqs overlap: Printing to file for hand checking')
                            with open('these_seqs_overlap_cluster.txt' , 'a') as a:
                                unresolved.append(list_of_GIs)
                                a.write(str(list_of_GIs) + '\n') 
                    else: 
            #get taxonomy for query(main species)
                        print("Parsing taxonomy for error sequences")
                        
                        ###here - is i tc_id or num_name?
                        hits = blast(i, list_of_GIs, c, taxdb)
                        #if theres only one lowest taxonomy hit, change
                        if hits.count(min(hits)) == 1:
                            finalseqs.add(str(two_dic[i][hit_levels.index(min(hits))]))
#                            print(str(two_dic[i][hit_levels.index(min(hits))]) + " had closer taxonomy hit")
                        else: #there are multiple lowest taxonomy hits
                            multfinalseqs.append(two_dic[i])
                            problem_dic[i] = two_dic[i]
 #                           print('Taxonomies had the multiple closest hits')
                else:
                    multfinalseqs.append(two_dic[i])
#                    print("Complement iden: " + str(iden) + " so pair is fine")
            else:
                multfinalseqs.append(two_dic[i])
#                print("Reverse Complement iden: " + str(iden) + " so pair is fine")
        else:
            multfinalseqs.append(two_dic[i])
#            print("High Aligned Identity: " + str(iden) + " so pair is fine")
    
            
    print("length of resolved = " + str(len(finalseqs)))
    print("length of choose multiple = " + str(len(multfinalseqs)))
    print("length of unresolved = " + str(len(unresolved)))
    
    with open("final_GIs.txt", "a") as o:
        for m in finalseqs:
            o.write(str(m)+"\n")
            
    with open("choose_mult.txt", "a") as o:
        for m in [num for pair in multfinalseqs for num in pair]:
            o.write(str(m)+"\n")

