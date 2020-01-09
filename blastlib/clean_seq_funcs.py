#!/usr/bin/env python
import time
from Bio import Entrez, SeqIO, AlignIO
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
from cleanlib.databasing import get_seqs_from_sqldb_GI, export_fasta, local_blast

import itertools
muscle_cline = MuscleCommandline(clwstrict=True)

def resolve_seqs(list_of_GIs, blastdb, gene, c):
    #give a list of GIs, checks the amount DNA/length, returns list of max
    amount_DNA = []
    length_DNA = []
    list_of_accs = []
    GI_to_pick = []
    iterator = get_seqs_from_sqldb_GI(list_of_GIs, "hseq", blastdb, gene, c)
    for seq_record in iterator:
        list_of_accs.append(seq_record.id)
        amount_DNA.append(seq_record.seq.count("A")+seq_record.seq.count("T")+seq_record.seq.count("C")+seq_record.seq.count("G"))
    for i in list_of_accs:
        for iter in c.execute("SELECT hit_length FROM blast WHERE accession = '"+ i +"';"):
            length_DNA.append(int(iter[0]))
    if amount_DNA.count(max(amount_DNA)) == 1:
        acc_to_pick = [list_of_accs[amount_DNA.index(max(amount_DNA))]]
    elif length_DNA.count(max(length_DNA)) == 1:
        acc_to_pick = [list_of_accs[length_DNA.index(max(length_DNA))]]
    else:
        try:
            sums = [amount_DNA[i] + length_DNA[i] for i in xrange(len(amount_DNA))] #python2
        except:
            sums = [amount_DNA[i] + length_DNA[i] for i in range(len(amount_DNA))] #python3
        acc_to_pick = [list_of_accs[i] for i, x in enumerate(sums) if x == max(sums)]
    for i in acc_to_pick:
        for iter in c.execute("SELECT GI from blast where accession = '"+i+"';"):
            GI_to_pick.append(str(iter[0]))
    return(GI_to_pick)

def alignment_reg(align_GIs, blastdb, qseqbool, gene, c):
    #qseqbool is a boolian that says whether or not the first GI is a qseq or not
    seqs = []
    if qseqbool == True:
        iterator = get_seqs_from_sqldb_GI(align_GIs[:1], "qseq", blastdb, gene, c)
        handle_string = StringIO()
        for seq in iterator:
            seqs.append(seq)
        iterator = get_seqs_from_sqldb_GI(align_GIs[1:], "hseq", blastdb, gene, c)
        handle_string = StringIO()
        for seq in iterator:
            seqs.append(seq)
    if qseqbool == False:
        iterator = get_seqs_from_sqldb_GI(align_GIs, "hseq", blastdb, gene, c)
        handle_string = StringIO()
        for seq in iterator:
            seqs.append(seq)
            
    SeqIO.write(seqs, handle_string, "fasta")
    data = handle_string.getvalue()
    stdout, stderr = muscle_cline(stdin=data)
    align = AlignIO.read(StringIO(stdout), "clustal")
    return(align)

def identity_calc(align):
    count = 0
    gaps = 0
    for col in range(0, len(align[0])):
        column = align[:,col].replace("N", "").replace("-", "")
        if len(column) > 1:
            if column[1:]==column[:-1]:
                count=count+1
        else:
            gaps=gaps+1
    if len(align[0])-gaps == 0:
        iden = 0.0
    else:
        iden = 100*(count/float((len(align[0])-gaps)))
    return(iden)

def alignment_rev_comp(align_GIs, blastdb, qseqbool, gene, c):
    #qseqbool is a boolian that says whether or not the first GI is a qseq or not
    seqs = []
    if qseqbool == True:
        iterator = get_seqs_from_sqldb_GI(align_GIs[:1], "qseq", blastdb, gene, c)
    if qseqbool == False:
        iterator = get_seqs_from_sqldb_GI(align_GIs[:1], "hseq", blastdb, gene, c)
    handle_string = StringIO()
    for seq in iterator:
        seqs.append(seq.reverse_complement())
    iterator = get_seqs_from_sqldb_GI(align_GIs[1:], "hseq", blastdb, gene, c)
    handle_string = StringIO()
    for seq in iterator:
        seqs.append(seq)
    SeqIO.write(seqs, handle_string, "fasta")
    data = handle_string.getvalue()
    stdout, stderr = muscle_cline(stdin=data)
    align = AlignIO.read(StringIO(stdout), "clustal")
    return(align)

def alignment_comp(align_GIs, blastdb, qseqbool, gene, c):
    #qseqbool is a boolian that says whether or not the first GI is a qseq or not
    seqs = []
    if qseqbool == True:
        iterator = get_seqs_from_sqldb_GI(align_GIs[:1], "qseq", blastdb, gene, c)
    if qseqbool == False:
        iterator = get_seqs_from_sqldb_GI(align_GIs[:1], "hseq", blastdb, gene, c)
    handle_string = StringIO()
    for seq in iterator:
        seq_comp = SeqRecord(seq.seq.complement())
        seq_comp.id = seq.id
        seq_comp.description = seq.description
        seqs.append(seq_comp)
    iterator = get_seqs_from_sqldb_GI(align_GIs[1:], "hseq", blastdb, gene, c)
    handle_string = StringIO()
    for seq in iterator:
        seqs.append(seq)
    SeqIO.write(seqs, handle_string, "fasta")
    data = handle_string.getvalue()
    stdout, stderr = muscle_cline(stdin=data)
    align = AlignIO.read(StringIO(stdout), "clustal")
    return(align)

def blast_all(blast_list, blast_nums, blast_tcids, c, gene, taxdb, blastdb):
    hit_levels_return = []
    iterator = get_seqs_from_sqldb_GI(blast_list, "hseq", blastdb, gene, c)
    export_fasta(iterator, "error_seqs.fa")
    local_blast(gene+"_db.fa", "error_seqs.fa")
    with open("error_seqs.fa.xml") as p:
        print("Parsing blast output")
        blast_records = NCBIXML.parse(p)
        #this iterates in a wierd way - note next at the end - this is due to the structure of ncbi's blast returns
        for rec in blast_records:
            #get GI - blast doesnt pull down GI nums anymore
            queryacc = str(rec.query.split()[0])
            for iter in c.execute("SELECT GI FROM blast WHERE accession = '" + queryacc + "'"):
                rec_GI = str(iter[0])
            #get taxonomy for tc_id
            rank = ''
            taxonomy = []
            sp_tc_id = blast_tcids[rec_GI]
            while rank != 'Family':
                for iter in c.execute("SELECT r.namestr, tc.parent_id FROM taxon_concepts tc, ranks r WHERE tc_id = '" + str(sp_tc_id) + "' AND tc.rank_id = r.rank_id"):
                    rank = iter[0]
                    sp_tc_id = iter[1]
                for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt WHERE ntt.name_id = n.name_id AND tc_id = '" + str(sp_tc_id) + "'"):
                    taxonomy.append(str(iter[0]))
            hit_levels_all = []
            for i in range(blast_nums[rec_GI]):
                query_acc = rec.query.split()[0]
    #            print(query_acc)
                hitdic = {}
                hitSp = []
                warning_counts = 0
                for alignment in rec.alignments:
                    for hsp in alignment.hsps:
                        identity=float(hsp.identities)/float(hsp.align_length)
                    if alignment.hit_def == query_acc:
                        pass
                    else:
                        for iter in c.execute("SELECT GI FROM blast WHERE accession='" + alignment.hit_def + "'"):
                            hitGI = (str(iter[0]))
                        hitdic[str(hitGI)] = identity
                #hitdic is GI: identity of all hits for a query
                if len(hitdic.values()) > 5:
                    maxiden = sorted(hitdic.values(), reverse = True)[0:5]
                else:
                    maxiden = hitdic.values()
                try:
                    hitGIs = [GI for GI, iden in hitdic.iteritems() if iden in maxiden] #python2
                except:
                    hitGIs = [GI for GI, iden in hitdic.items() if iden in maxiden] #python3
                #hitGIs - GIs with the highest identity
                for hitGI in hitGIs:
                    prev_len = len(hitSp)
        #get species for all the top hits for each GI - hitsp
                    for iter in c.execute("SELECT Species FROM blast WHERE GI='" + hitGI + "'"):
                        hitSp.append(str(iter[0]))
                    if len(hitSp) == prev_len:
                        warning_counts += 1
                        #print("Warning: " + str(hitGI) + " is a top hit, but is not in the blast database")
                hit_levels = []
                if len(hitSp) > 0:
                    if warning_counts < 2:
                        for hit_species in hitSp:
                            hit = -1
                            for iter in c.execute("SELECT ntt.tc_id FROM names_to_taxonconcepts ntt, names n WHERE n.name_id = ntt.name_id AND namestr='" + hit_species + "'"):
                                sp_tc_id = iter[0]
                            while hit == -1:
                                for iter in c.execute("SELECT r.namestr, tc.parent_id, tc.tc_id FROM taxon_concepts tc, ranks r WHERE tc_id = '" + str(sp_tc_id) + "'AND tc.rank_id = r.rank_id"):
                                    sp_tc_id = iter[1]
                                for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt WHERE ntt.name_id = n.name_id AND tc_id = '" + str(sp_tc_id) + "'"):
                                    if str(iter[0]) in taxonomy:
                                        hit = taxonomy.index(iter[0])
                                        hit_levels.append(hit)
                        hit_levels_all.append(min(hit_levels))
                    else:
    #                    print('Warning: GI has two or more hits that aren\'t in the database, will ignore')
                        hit_levels_all.append(1000)
                else:
    #                print('Warning: no taxonomy info for current GI\'s hits - will look at other possible GI')
                    hit_levels_all.append(1000)
                if i != blast_nums[rec_GI]-1:
                    rec = next(blast_records)
            hit_levels_return.append(hit_levels_all)
    return(hit_levels_return)
    
    
def blast(sp_tc_id, align_GIs, c, taxdb):
    rank = ''
    taxonomy = []
    ent_query = get_blast_query(taxdb)
    while rank != 'Family':
        for iter in c.execute("SELECT r.namestr, tc.parent_id FROM taxon_concepts tc, ranks r WHERE tc_id = '" + str(sp_tc_id) + "'AND tc.rank_id = r.rank_id"):
            rank = iter[0]
            sp_tc_id = iter[1]
        for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt WHERE ntt.name_id = n.name_id AND tc_id = '" + str(sp_tc_id) + "'"):
            taxonomy.append(str(iter[0]))
    hit_levels_all = []
# do blast of the contrasting GIs
#       have to repull down because data is removed from handle and seqs and can't parse back from a string
    for query_GI in align_GIs:
#        print("Blasting " + query_GI + " to find the closest hit taxonomy")
        error = True
        while error == True:
            try:
                result_handle = NCBIWWW.qblast("blastn", "nt", query_GI, entrez_query=ent_query, word_size=28, hitlist_size=100)
                error = False
            except:
                print("Error, trying again")
                time.sleep(10)
        blast_records = NCBIXML.read(result_handle)
        hitdic = {}
        hitSp = []
        warning_counts = 0
        for alignment in blast_records.alignments:
            for hsp in alignment.hsps:
                identity=float(hsp.identities)/float(hsp.align_length)
            if alignment.title.split("|")[1] == query_GI:
                pass
            else:
                hitdic[str(alignment.title.split("|")[1])] = identity
        if len(hitdic.values()) > 5:
            maxiden = sorted(hitdic.values(), reverse = True)[0:5]
        else:
            maxiden = hitdic.values()
        try:
            hitGIs = [GI for GI, iden in hitdic.iteritems() if iden in maxiden] #python2
        except:
            hitGIs = [GI for GI, iden in hitdic.items() if iden in maxiden] #python3
        for hitGI in hitGIs:
            prev_len = len(hitSp)
#get species for all the top hits for each GI
            for iter in c.execute("SELECT Species FROM blast WHERE GI='" + hitGI + "'"):
                hitSp.append(str(iter[0]))
            if len(hitSp) == prev_len:
                warning_counts += 1
#                print("Warning: " + str(hitGI) + " is a top hit, but is not in the blast database")
        hit_levels = []
# for each hit species, find where taxonomy overlaps with query, want the GI of the lowest
        if len(hitSp) > 0:
            if warning_counts < 2:
                for hit_species in hitSp:
                    hit = -1
                    for iter in c.execute("SELECT ntt.tc_id FROM names_to_taxonconcepts ntt, names n WHERE n.name_id = ntt.name_id AND namestr='" + hit_species + "'"):
                        sp_tc_id = iter[0]
                    while hit == -1:
                        for iter in c.execute("SELECT r.namestr, tc.parent_id, tc.tc_id FROM taxon_concepts tc, ranks r WHERE tc_id = '" + str(sp_tc_id) + "'AND tc.rank_id = r.rank_id"):
                            sp_tc_id = iter[1]
                        for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt WHERE ntt.name_id = n.name_id AND tc_id = '" + str(sp_tc_id) + "'"):
                            if str(iter[0]) in taxonomy:
                                hit = taxonomy.index(iter[0])
                                hit_levels.append(hit)
                                print(hit)
                print(hit_levels)
                hit_levels_all.append(min(hit_levels))
            else:
#                print('Warning: GI has two or more hits that aren\'t in the database, will ignore')
                hit_levels_all.append(1000)
        else:
#            print('Warning: no taxonomy info for ' + str(query_GI) + '\'s hits - will choose other possible GI')
            hit_levels_all.append(1000)
    return(hit_levels_all)

def tiling(list_of_GIs_local, gene, blastdb, c):
    idens_local = []
    start_stop_local = []
    for m in list_of_GIs_local:
        for iter in c.execute("SELECT GI from blast WHERE Gene_name = '"+gene+"' LIMIT 1;"):
            qseqGI = str(iter[0])
        GIs_to_align = [qseqGI, m]
        #need to change this
        idens_for_ind = []
        alignment = alignment_reg(GIs_to_align, blastdb, True, gene, c)
        iden = identity_calc(alignment)
        idens_for_ind.append(iden)
        if iden < 70:
            alignment = alignment_rev_comp(GIs_to_align, blastdb, True, gene, c)
            iden = identity_calc(alignment)
            idens_for_ind.append(iden)
            if iden < 70: 
                alignment = alignment_comp(GIs_to_align, blastdb, True, gene, c)
                iden = identity_calc(alignment)
                idens_for_ind.append(iden)
        span = 0
        #get start
        for i in range(len(alignment[0])):
            col = alignment[:, i]
            if '-' not in col:
                span += 1
            if span == 10:
                break
            elif span > 0 and '-' in col:
                span = 0
        start = i-8
        span = 0
        #get stop
        for i in reversed(range(len(alignment[0]))):
            col = alignment[:, i]
            if '-' not in col:
                span += 1
            if span == 10:
                break
            elif span > 0 and '-' in col:
                span = 0
        end = i+10
        start_stop_local.append([start, end])
        idens_local.append(max(idens_for_ind))
    return(idens_local, start_stop_local)

def get_blast_query(taxonomy):
    
    import sqlite3
    conn = sqlite3.connect(taxonomy)
    c = conn.cursor()
    count = 100
    family_ranks = []
    #right now this just gets the families, but can be edited so can choose ranks
    #while count > 50:
    for iter in c.execute("SELECT n.namestr FROM names n, names_to_taxonconcepts ntt, taxon_concepts tc, ranks r WHERE ntt.name_id = n.name_id AND ntt.tc_id = tc.tc_id AND tc.rank_id = r.rank_id AND r.namestr = 'Family'"):
        family_ranks.append(str(iter[0]))
        #count = len(family_ranks)
        #print(count)
    family_orgs = ["(" + i + "[Organism])" for i in family_ranks]
    query = " OR ".join(family_orgs)
    return(query)


