#!/usr/bin/env python
import time
from Bio import Entrez, SeqIO, AlignIO
from StringIO import StringIO
from Bio.Align.Applications import MuscleCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML
import itertools
muscle_cline = MuscleCommandline(clwstrict=True)
Entrez.email = "sunray1@ufl.edu"

def resolve_seqs(list_of_GIs):
    #give a list of GIs, checks the amount DNA/length, returns list of max
    joined_GIs = ",".join(list_of_GIs)
    amount_DNA = []
    length_DNA = []
    error = True
    while error == True:
        try:
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=joined_GIs)
            error = False
        except:
            print("Error, trying again")
            time.sleep(10)
    for seq_record in SeqIO.parse(handle, "fasta"):
        amount_DNA.append(seq_record.seq.count("A")+seq_record.seq.count("T")+seq_record.seq.count("C")+seq_record.seq.count("G"))
        length_DNA.append(len(seq_record))
#    print(amount_DNA)
#    print(length_DNA)
    if amount_DNA.count(max(amount_DNA)) == 1:
        GI_to_pick = [list_of_GIs[amount_DNA.index(max(amount_DNA))]]
    elif length_DNA.count(max(length_DNA)) == 1:
        GI_to_pick = [list_of_GIs[length_DNA.index(max(length_DNA))]]
    else:
        sums = [amount_DNA[i] + length_DNA[i] for i in xrange(len(amount_DNA))]
        GI_to_pick = [list_of_GIs[i] for i, x in enumerate(sums) if x == max(sums)]
    return(GI_to_pick)

def alignment_reg(align_GIs):
    joined_GIs = ",".join(align_GIs)
    error = True
    while error == True:
        try:
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=joined_GIs)
            error = False
        except:
            print("Error, trying again")
            time.sleep(10)
    seqs =  SeqIO.parse(handle, "fasta")
    handle_string = StringIO()
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

def alignment_rev_comp(align_GIs):
    firstGI = align_GIs[0]
    error = True
    while error == True:
        try:
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=firstGI)
            error = False
        except:
            print("Error, trying again")
            time.sleep(10)
    seqs =  SeqIO.parse(handle, "fasta")
    seq = next(seqs).reverse_complement()
    handle_string = StringIO()
    SeqIO.write(seq, handle_string, "fasta")
    otherGIs = align_GIs[1:]
    error = True
    while error == True:
        try:
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=otherGIs)
            error = False
        except:
            print("Error, trying again")
            time.sleep(10)
    seqs =  SeqIO.parse(handle, "fasta")
    SeqIO.write(seqs, handle_string, "fasta")
    data = handle_string.getvalue()
    stdout, stderr = muscle_cline(stdin=data)
    align = AlignIO.read(StringIO(stdout), "clustal")
    return(align)

def alignment_comp(align_GIs):
    firstGI = align_GIs[0]
    error = True
    while error == True:
        try:
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=firstGI)
            error = False
        except:
            print("Error, trying again")
            time.sleep(10)
    seqs =  SeqIO.parse(handle, "fasta")
    seq = next(seqs)
    #since seqrecord doesnt have a complement function yet, we have to call the complement on the Seq object(the sequence itself) and create a seqrecord
    seq_comp = SeqRecord(seq.seq.complement())
    seq_comp.id = seq.id
    seq_comp.description = seq.description
    handle_string = StringIO()
    SeqIO.write(seq_comp, handle_string, "fasta")
    otherGIs = align_GIs[1:]
    error = True
    while error == True:
        try:
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=otherGIs)
            error = False
        except:
            print("Error, trying again")
            time.sleep(10)
    seqs =  SeqIO.parse(handle, "fasta")
    SeqIO.write(seqs, handle_string, "fasta")
    data = handle_string.getvalue()
    stdout, stderr = muscle_cline(stdin=data)
    align = AlignIO.read(StringIO(stdout), "clustal")
    return(align)

def blast_all(blast_list, blast_nums, blast_tcids, c):
    hit_levels_return = []
    joined_GIs = ",".join(blast_list)
    error = True
    while error == True:
        try:
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=joined_GIs)
            error = False
        except:
            print("Error, trying again")
            time.sleep(10)
    seqs =  SeqIO.parse(handle, "fasta")
    handle_string = StringIO()
    SeqIO.write(seqs, handle_string, "fasta")
    data = handle_string.getvalue()
    error = True
    while error == True:
        try:
            result_handle = NCBIWWW.qblast("blastn", "nt", data, entrez_query='((Papilionoidea[Organism]) OR Hedylidae[Organism]) OR Hesperiidae[Organism]', word_size=28, hitlist_size=100)
            error = False
        except:
            print("Error, trying again")
            time.sleep(10)
    blast_records = NCBIXML.parse(result_handle)
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
                if alignment.title.split("|")[3] == query_acc:
                    pass
                else:
                    hitdic[str(alignment.title.split("|")[1])] = identity
            if len(hitdic.values()) > 5:
                maxiden = sorted(hitdic.values(), reverse = True)[0:5]
            else:
                maxiden = hitdic.values()
            hitGIs = [GI for GI, iden in hitdic.iteritems() if iden in maxiden]
            for hitGI in hitGIs:
                prev_len = len(hitSp)
    #get species for all the top hits for each GI
                for iter in c.execute("SELECT Species FROM blast WHERE GI='" + hitGI + "'"):
                    hitSp.append(str(iter[0]))
                if len(hitSp) == prev_len:
                    warning_counts += 1
#                    print("Warning: " + str(hitGI) + " is a top hit, but is not in the blast database")
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
    
    
def blast(sp_tc_id, align_GIs, c):
    rank = ''
    taxonomy = []
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
                result_handle = NCBIWWW.qblast("blastn", "nt", query_GI, entrez_query='((Papilionoidea[Organism]) OR Hedylidae[Organism]) OR Hesperiidae[Organism]', word_size=28, hitlist_size=100)
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
        hitGIs = [GI for GI, iden in hitdic.iteritems() if iden in maxiden]
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

def tiling(list_of_GIs_local, gene):
    idens_local = []
    start_stop_local = []
    for m in list_of_GIs_local:
        if gene == 'CAD':
            GIs_to_align = [m, '440201778']
        if gene == 'EF-1a':
            GIs_to_align = [m, '752854557']
        if gene == 'RpS5':
            GIs_to_align = [m, '342356392']
        if gene == 'ArgKin':
            GIs_to_align = [m, '161088179']
        if gene == 'CAT':
            GIs_to_align = [m, '389612591']
        if gene == 'CoA':
            GIs_to_align = [m, '345096916']
        if gene == 'DDC':
            GIs_to_align = [m, '315493443']
        if gene == 'GAPDH':
            GIs_to_align = [m, '914615339']
        if gene == 'HCL':
            GIs_to_align = [m, '298398979']
        if gene == 'IDH':
            GIs_to_align = [m, '914615459']
        if gene == 'MDH':
            GIs_to_align = [m, '914614510']
        if gene == 'RpS2':
            GIs_to_align = [m, '342356316']
        if gene == 'Wgl':
            GIs_to_align = [m, '58294437']
        if gene == 'COI_trnL_COII':
            GIs_to_align = [m, 'GU365907']
        #need to change this
        idens_for_ind = []
        alignment = alignment_reg(GIs_to_align)
        iden = identity_calc(alignment)
        idens_for_ind.append(iden)
        if iden < 70:
            alignment = alignment_rev_comp(GIs_to_align)
            iden = identity_calc(alignment)
            idens_for_ind.append(iden)
            if iden < 70: 
                alignment = alignment_comp(GIs_to_align)
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
        
