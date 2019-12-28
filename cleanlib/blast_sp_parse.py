#!/usr/bin/env python
#Script does self blast of those with multiple species/gene hits to make sure the top hit is the rejected hit
from blastlib.clean_seq_funcs import get_blast_query

def test_resolved_seqs(infile, blastdb, taxdb, email):
    ent_query = get_blast_query(taxdb)
    import sqlite3, sys, time
    from Bio.Blast import NCBIWWW, NCBIXML
    from blastlib.clean_seq_funcs import resolve_seqs, alignment_comp, alignment_reg, alignment_rev_comp, blast, blast_all, identity_calc, tiling
    from cleanlib.databasing import get_seqs_from_sqldb, export_fasta, create_blast_db, local_blast
    conn = sqlite3.connect(blastdb)
    c = conn.cursor()
    c.execute("ATTACH '" + taxdb + "' as 'tax'")
    error_dic = {}
    blast_dic_nums = {}
    blast_dic_tcids = {}
    seqs_to_blast = []
    finalseqs = set()
    multseqs = []
    
    
    #do self blast and print to file
    with open(infile) as fasta_file:
        print("Blasting " + infile)
        records = fasta_file.read()
    numqueries = records.count('>')
    # 
    # #get all seqs that have that gene to make db out of
    gene = infile.split("_accession_nums")[0]
    access_list = set()
    for iter in c.execute("SELECT accession FROM blast WHERE Gene_name = '"+ gene +"';"):
        access_list.add(iter[0])
    #print(len(access_list))
    iterator = get_seqs_from_sqldb(access_list, "hseq", blastdb)
    export_fasta(iterator, gene+"_db.fa")
    create_blast_db(gene+"_db.fa")
    local_blast(gene+"_db.fa", infile)
    
    #open self blast file for parsing
    #make dictionary of query species: query GI of those that don't have the top hit as the same species
    
    with open(infile+".xml") as p:
        print("Parsing blast output")
        blast_recs=NCBIXML.parse(p)
        count = 0
        for rec in blast_recs:
            count += 1
            print(str(round(float(count)/float(numqueries)*100, 2))+ "%")
#figure out a new way to do this
            queryAcc = str(rec.query.split()[0])
            for iter in c.execute("SELECT GI FROM blast WHERE accession='" + queryAcc + "';"):
                queryGI = (str(iter[0]))
            #hitdic is GIs and idens of all in .xml
            hitdic = {}
            hitSp = set()
            for alignment in rec.alignments:
                for hsp in alignment.hsps:
                    identity=float(hsp.identities)/float(hsp.align_length)
                    #print(identity)
                if alignment.title.split("|")[1] == queryGI:
                    pass
                else:
                    hitdic[str(alignment.title.split("|")[1])] = identity
            #print(hitdic)
            maxiden = max(hitdic.values())
            try:
                hitGIs = [GI for GI, iden in hitdic.iteritems() if iden == maxiden] #python2
            except:
                hitGIs = [GI for GI, iden in hitdic.items() if iden == maxiden] #python3
            for iter in c.execute("SELECT tc_id FROM blast WHERE GI='" + queryGI + "'"):
                querySp = (str(iter[0]))
            for i in hitGIs:
                for iter in c.execute("SELECT tc_id FROM blast WHERE GI='" + i + "'"):
                    hitSp.add(str(iter[0]))
            #only look at top 5 if it doesnt hit the top hit - faster
            #sort hitdic by idens
            if querySp not in hitSp:
                hitSp = set()
                if len(hitdic.values()) > 5:
                    maxiden = sorted(hitdic.values(), reverse = True)[0:5]
                else:
                    maxiden = hitdic.values()
                try:
                    hitGIs = [GI for GI, iden in hitdic.iteritems() if iden in maxiden] #python2
                except:
                    hitGIs = [GI for GI, iden in hitdic.items() if iden in maxiden] #python3
                #get hit species (actually tc_ids)
                for i in hitGIs:
                    for iter in c.execute("SELECT tc_id FROM blast WHERE GI='" + i + "'"):
                        hitSp.add(str(iter[0]))
                #if species of query not in the list of hit species
                if querySp not in hitSp:
                    error_dic[querySp] = queryGI
                else:
                    c.execute("UPDATE blast SET decision='Longest or most info, good top hit/chosen' WHERE GI='" + queryGI + "';")
                    finalseqs.add(queryGI)
            else:
                c.execute("UPDATE blast SET decision='Longest or most info, good top hit/chosen' WHERE GI='" + queryGI + "';")
                finalseqs.add(queryGI)
    count = 0
    #error_dic is dictionary that has species = GI that doesn't match species when self blast
    #error_dic['21204'] = '316994286'
    ##go through error dictionary and align the 'same' gene/species to see if they're weird looking
    print("Checking nonmatching sequences")
    for tc_id in error_dic:
        count += 1
        print(str(round(float(count)/float(len(error_dic))*100, 2))+"%")
        list_of_GIs = []
        for iter in c.execute("SELECT Gene_name FROM blast WHERE GI ='" + error_dic[tc_id] + "'"):
            gene_name = str(iter[0])
        for iter in c.execute("SELECT GI FROM blast WHERE tc_id = '" + tc_id + "' and Gene_name = '" + gene_name + "'"):                
            list_of_GIs.append(str(iter[0]))
        first_GI = list_of_GIs[0]
        other_GIs = list_of_GIs[1:]
        pair_check = []
        #align each seq to the first one - first one acts as the 'default' direction
        #print(other_GIs)
        for x, GI in enumerate(other_GIs):
            GI_pair = [first_GI, GI]
            alignment = alignment_reg(GI_pair, blastdb, False)
            iden = identity_calc(alignment)
            if iden < 90:
#                print("Low Aligned Identity: " + str(iden))
                alignment = alignment_rev_comp(GI_pair, blastdb, False)
                iden = identity_calc(alignment)
                if iden < 90: 
#                    print("Low Reverse Complement Aligned Identity: " + str(iden))
                    alignment = alignment_comp(GI_pair, blastdb, False)
                    iden = identity_calc(alignment)
                    if iden < 90:
#                        print("Low Complement Aligned Identity: " + str(iden))
                        pair_check.append(0)
                    else:
                        pair_check.append(1)
#                        print("Complement iden: " + str(iden) + " so pair is fine")
                else:
                    pair_check.append(1)
#                    print("Reverse Complement iden: " + str(iden) + " so pair is fine")
            else:
                pair_check.append(1)
#                print("High Aligned Identity: " + str(iden) + " so pair is fine")
        #print(pair_check)
        #1 is when the GI aligned to chosen one just fine, 0 is when it doesn't
        if all(i == 1 for i in pair_check):
            c.execute("UPDATE blast SET decision='Sequence did not have same top blast species, but all aligned correctly/Chosen' WHERE GI='" + error_dic[tc_id] + "';")
            finalseqs.add(error_dic[tc_id])
        else:
            idens, start_stop = tiling(list_of_GIs, gene_name, blastdb, c)
            current_start = -1
            current_stop = -1 
            result = []
            if all(i > 70 for i in idens):
                for start, stop in sorted(start_stop):
                    if start > current_stop:
                        result.append((start, stop))
                        current_start, current_stop = start, stop
                    else:
                        current_stop = max(current_stop, stop)
                        result[-1] = (current_start, current_stop)
                if len(result) == len(start_stop):
#                    print("Seqs align to different regions of probe, choosing all")
                    for x in list_of_GIs:
                        c.execute("UPDATE blast SET decision='Best sequence in tiling analysis/Chosen' WHERE GI='" + x + "';")
                        finalseqs.add(x)
                else:
#                    print('Seqs overlap: Printing to file for hand checking')
                    #these sequences are wierd because they align to the whole, but don't align together well
                    with open('these_seqs_overlap.txt' , 'a') as a:
                        a.write(str(list_of_GIs) + '\n') 
            else:
#                print('Somethings up with a sequence - printing to check - will blast')
                pair_check.append(0)
                blast_dic_nums[list_of_GIs[0]] = len(list_of_GIs)
                blast_dic_tcids[list_of_GIs[0]] = tc_id
                seqs_to_blast.append(list_of_GIs)
                with open('seqs_to_be_blasted.txt', 'a') as a:
                    a.write(str(list_of_GIs) + '\n')

    if len(seqs_to_blast) > 0:
        print("Blasting error sequences (seqs do not align together and one doesn't align to whole)")        
        seqs_to_blast_flat = [item for sublist in seqs_to_blast for item in sublist]
        try:
            hits_all = blast_all(seqs_to_blast_flat, blast_dic_nums, blast_dic_tcids, c, email, taxdb)
        except:
            time.sleep(5)
            hits_all = blast_all(seqs_to_blast_flat, blast_dic_nums, blast_dic_tcids, c, email, taxdb)
#        print(hits_all)
        print("Parsing taxonomy for error sequences")
        for x, list_of_GIs in enumerate(seqs_to_blast):
            hits = hits_all[x]
            tc_id = blast_dic_tcids[list_of_GIs[0]]
        #if theres only one lowest taxonomy hit and its not itself, change
            if hits.count(min(hits)) == 1 and error_dic[tc_id] != list_of_GIs[hits.index(min(hits))]:
                finalseqs.add(list_of_GIs[hits.index(min(hits))])
#                print(str(list_of_GIs[hits.index(min(hits))]) + " had closer taxonomy hit")
            elif hits.count(min(hits)) == 1 and error_dic[tc_id] == list_of_GIs[hits.index(min(hits))]:
                finalseqs.add(error_dic[tc_id])
#                print(str(list_of_GIs[hits.index(min(hits))]) + " was previously chosen")
            else: #there are multiple lowest taxonomy hits
#                print('Taxonomies had multiple closest hits')
                index_pos = []
                count = 0
                for x in hits:
                    if x == min(hits):
                        index_pos.append(count)
                    count+=1
                mult_GIs = [list_of_GIs[x] for x in index_pos]
                GI_to_pick = resolve_seqs(mult_GIs)
                #if theres only one chosen and it wasnt the one already picked...add to change dic
                if len(GI_to_pick) == 1 and error_dic[tc_id]!=GI_to_pick[0]:
                    finalseqs.add(GI_to_pick[0])
#                    print(str(GI_to_pick) + " chosen")
                #if theres only one max length and it was already picked
                elif len(GI_to_pick) == 1 and error_dic[tc_id]==GI_to_pick[0]:
#                    print(str(GI_to_pick) + " was previously chosen")
                    finalseqs.add(GI_to_pick[0])
                else:
                #this only happens if the originally picked one is crappy and the rest are the same
#                    print("Multiple choices: " + str(GI_to_pick))
                    multseqs.append(error_dic[tc_id])
                    #Go to cluster analysis
            
    print('length of resolved=' + str(len(finalseqs)))
    print('length of not resolved = ' + str(len(multseqs)))
            
    with open("final_GIs.txt", "a") as o:
       for m in finalseqs:
           o.write(str(m)+"\n")
    
    #to cluster
    with open("multiple_gene_choices.txt", "a") as o:
        for m in multseqs:
            o.write(str(m)+"\n")
        
                    
                    
    conn.commit()
    conn.close()
