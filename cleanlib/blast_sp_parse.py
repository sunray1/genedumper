#!/usr/bin/env python
#Script does self blast of those with multiple species/gene hits to make sure the top hit is the rejected hit
from blastlib.clean_seq_funcs import get_blast_query


def first_blast(blastdb):
    #do self blast and print to file
    import sqlite3
    from cleanlib.databasing import get_seqs_from_sqldb, export_fasta, create_blast_db, local_blast
    genes = set()
    conn = sqlite3.connect(blastdb)
    c = conn.cursor()
    # #get all seqs that have that gene to make db out of
    for iter in c.execute("SELECT Gene_name from blast GROUP BY Gene_name;"):
        genes.add(iter[0])
    for gene in genes:
        #make query
        print("Blasting " + gene +" sequences")
        access_list_qseq = set()
        for iter in c.execute("SELECT accession FROM blast WHERE decision = 'Longest or most info (not checked)/Chosen' and Gene_name = '"+gene+"';"):
            access_list_qseq.add(str(iter[0]))
        for iter in c.execute("SELECT accession FROM blast WHERE decision = 'Only or best choice in tiling analysis/chosen'and Gene_name = '"+gene+"';"):
            access_list_qseq.add(str(iter[0]))
        iterator = get_seqs_from_sqldb(access_list_qseq, "hseq", blastdb, c)
        export_fasta(iterator, gene+"_qseqs.fa")
        #make db and blast
        if len(access_list_qseq) != 0:
            access_list = set()
            for iter in c.execute("SELECT accession FROM blast WHERE Gene_name = '"+ gene +"' AND tc_id != '0';"):
                access_list.add(str(iter[0]))
            #print(len(access_list))
            iterator = get_seqs_from_sqldb(access_list, "hseq", blastdb, c)
            export_fasta(iterator, gene+"_db.fa")
            create_blast_db(gene+"_db.fa")
            local_blast(gene+"_db.fa", gene+"_qseqs.fa")
        
def test_resolved_seqs(infile, blastdb, taxdb, email):
    ent_query = get_blast_query(taxdb)
    import sqlite3, sys, time, subprocess
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
    #open self blast file for parsing
    #make dictionary of query species: query GI of those that don't have the top hit as the same species
    
    with open(infile) as p:
        print("Parsing reciprocal blast output from " + infile)
        blast_recs=NCBIXML.parse(p)
        count = 0
        numqueries = subprocess.check_output("grep -c '>' "+ infile.split(".xml")[0], shell=True)
        for rec in blast_recs:
            count += 1
            print(str(round(float(count)/float(numqueries)*100, 2))+ "%")
#figure out a new way to do this
            queryAcc = str(rec.query.split()[0])
            #this broke when ncbi updated
            for iter in c.execute("SELECT GI FROM blast WHERE accession='" + queryAcc + "';"):
                queryGI = (str(iter[0]))
            #print(queryGI)
            #hitdic is GIs and idens of 20 in .xml for each rec in blast_rec
            hitdic = {}
            hitSp = set()
            count1 = 0
            for alignment in rec.alignments:
                    for hsp in alignment.hsps:
                        identity=float(hsp.identities)/float(hsp.align_length)
                        #print(identity)
                    if alignment.hit_def == queryAcc:
                        pass
                    else:
                        for iter in c.execute("SELECT GI FROM blast WHERE accession='" + alignment.hit_def + "'"):
                            hitGI = (str(iter[0]))
                        count1 += 1
                        hitdic[str(hitGI)] = identity
                        if count1 >= 20:
                            break
                # if alignment.title.split("|")[1] == queryGI:
                #     pass
                # else:
                #     hitdic[str(alignment.title.split("|")[1])] = identity
            
            if len(hitdic.values()) == 0:
                #####################chose next one#######################
                print('No matching hits from blast')
                c.execute("UPDATE blast SET Decision='Chosen, but then did not blast to anything/Not chosen' WHERE GI = '" + queryGI + "';")
                pass
            else:
                ##first looks at just the top iden value
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
                ##if it doesn't work, look at top 5 idens - ignore top hit b/c already looked at
                if querySp not in hitSp:
                    hitSp = set()
                    if len(hitdic.values()) > 5:
                        maxiden = sorted(set(hitdic.values()), reverse = True)[1:5]
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
#################will rewrite over if multiple genes per querysp######################
####################can't align all together if they align to different parts####################
                    if querySp not in hitSp:
                        c.execute("UPDATE blast SET decision='Chosen, but does not reciprocal blast' WHERE GI='" + queryGI + "';")
                        error_dic[querySp] = [queryGI]
                        
                    else:
                        c.execute("UPDATE blast SET decision='Longest or most info, good top hit/chosen' WHERE GI='" + queryGI + "';")
                else:
                    c.execute("UPDATE blast SET decision='Longest or most info, good top hit/chosen' WHERE GI='" + queryGI + "';")
    count = 0
    #error_dic is dictionary that has species = GI that doesn't match species when self blast
    # error_dic = {'81': '908332352', '56': '227435787', '60': '224939202', '39': '497153612', '157': '224939198', '87': '695134008', '628': '695134046', '184': '1174533403', '116': '558477133', '121': '316994106', '206': '545690375', '429': '1160414548', '431': '227436029', '480': '695134034', '105': '545690685', '651': '695134010', '413': '529377043', '373': '443611395', '252': '529388651', '354': '529386663', '335': '443611449', '584': '295078734', '375': '1559462858', '523': '1051545064', '350': '443611513', '410': '302487953', '527': '316994100', '244': '575502624', '355': '443611465', '595': '317467717', '380': '443611323', '357': '529374537', '577': '317467813', '363': '675401803', '406': '443611367', '372': '73533766', '656': '299833006'}
    ##go through error dictionary and align the 'same' gene/species to see if they're weird looking
    print("Aligning sequences where chosen does not reciprocal blast to another")
    for tc_id in error_dic:
        count += 1
        print(str(round(float(count)/float(len(error_dic))*100, 2))+"%")
        list_of_GIs = []
        for iter in c.execute("SELECT Gene_name FROM blast WHERE GI ='" + error_dic[tc_id] + "'"):
            gene = str(iter[0])
        for iter in c.execute("SELECT GI FROM blast WHERE tc_id = '" + tc_id + "' and Gene_name = '" + gene + "'"):                
            list_of_GIs.append(str(iter[0]))
        first_GI = list_of_GIs[0]
        other_GIs = list_of_GIs[1:]
        pair_check = []
        #align each seq to the first one - first one acts as the 'default' direction
        for x, GI in enumerate(other_GIs):
            GI_pair = [first_GI, GI]
            alignment = alignment_reg(GI_pair, blastdb, False, gene, c)
            iden = identity_calc(alignment)
            if iden < 90:
#                print("Low Aligned Identity: " + str(iden))
                alignment = alignment_rev_comp(GI_pair, blastdb, False, gene, c)
                iden = identity_calc(alignment)
                if iden < 90: 
#                    print("Low Reverse Complement Aligned Identity: " + str(iden))
                    alignment = alignment_comp(GI_pair, blastdb, False, gene, c)
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
        else:
            ##either all are 0s if first one is wrong or one is zero where one is wrong
#           print('Somethings up with a sequence - printing to check - will blast')
            seqs_to_blast.append(list_of_GIs)
            list_of_GIs_str = str(list_of_GIs).replace("[", "(").replace("]", ")")
            c.execute("UPDATE blast SET Decision='Chosen, but does not reciprocal blast and all do not align' WHERE GI IN " + list_of_GIs_str + ";") 

    # seqs_to_blast = [['695134046', '317467897'], ['316994106', '295069424', '414303142'], ['695134034', '71836048'], ['302487953', '73533700', '66096804'], ['575502624', '1723904125', '1723904119'], ['73533766', '443611441', '443611439'], ['296792260', '1531247074', '1531246852', '1531245174', '1531244766', '1531244542', '1531244412', '870902494', '520760133', '520760131', '331691381', '331690706', '331690580', '331690574', '331690568', '296467088', '296467084', '331691379', '331690792', '974999640', '974999638', '1531247630', '1531247410', '1531247308', '1531247180', '1531247130', '1531246686', '1531246650', '1531246450', '1531246148', '1531246132', '1531245926', '1531245702', '1531245624', '1531245592', '1531245462', '1531245454', '1531245440', '1531245354', '1531244888', '1531244834', '1531244804', '1531244664', '1531244614', '1531244586', '1531244178', '1531244132', '1531244016', '520760163', '520760159', '520760149', '331691355', '331690624', '331690576', '331690564', '331690544', '307641798', '296469132', '296468960', '296467194', '296466760', '304270861', '331690704', '63030162', '974999636', '227436371', '1532637164', '1532637082', '1532637024', '1531247684', '1531247306', '1531247040', '1531246466', '1531246208', '1531246108', '1531245848', '1531245842', '1531245638', '1531245614', '1531245574', '1531245376', '1531244734', '1531244658', '1531244080', '1531244012', '1531243910', '633896164', '520760155', '520759941', '331691337', '331690536', '300203102', '156619418', '156619410', '974999630', '296467200', '304270883', '1532637150', '1532637144', '1532637064', '1532637028', '1532637018', '730837173', '633896166', '520760157', '520760151', '331691341', '331691339', '331690582', '331690562', '1557923373', '1557923370', '156619416', '156619412', '1532637074', '1532637022', '1532637016', '156619414', '156619408', '730837157', '633896160', '1532637100', '1532637090', '1532637072', '1532637038', '633896162', '1042273875', '520760209', '313173207', '296468852', '88604889', '1532637050', '520760225', '520760219', '520760213', '296468054', '1532637044', '520760221', '520760193', '296792296', '312928999', '1557923367', '1557923364', '156619420', '331690820', '156619406', '520760211', '227436369', '296463378', '296463050', '227436365', '227436367', '296463374', '870901906', '870902242', '299833006', '76008949', '296464328', '296465368', '296464998', '296465182', '296464956', '296465198', '296463914']]
    ##sequences don't align together and chosen hit doesn't blast to same species (one area of tiling)
    if len(seqs_to_blast) > 0:
        for iter in c.execute("SELECT Gene_name FROM blast WHERE GI ='" + seqs_to_blast[0][0] + "'"):
            gene = str(iter[0])
        seqs_to_blast_flat = [item for sublist in seqs_to_blast for item in sublist]
        #queryGI: nearest tax rank of hit
        nearest_hit_taxo_dic = blast_all(seqs_to_blast_flat, c, gene, taxdb, blastdb)
        print("Parsing taxonomy for error sequences")
        for i in seqs_to_blast:
            hit_taxonomies = [nearest_hit_taxo_dic[x] for x in i]
            #want max tcid - that is closest to query
            maxtcid = max(hit_taxonomies)
            GI_to_pick = [GI for GI in enumerate(hit_taxonomies) if iden == maxtcid]
            if len(GI_to_pick) == 1:
                c.execute("UPDATE blast SET decision='Blast hits have closest taxonomy/Chosen' WHERE GI='" + GI_to_pick[0] + "';")
            else:
                GI_to_pick_str = str(GI_to_pick).replace("[", "(").replace("]", ")")
                c.execute("UPDATE blast SET Decision='To cluster analysis/Chosen' WHERE GI IN " + GI_to_pick_str + ";")    
                    
                    
    conn.commit()
    conn.close()
