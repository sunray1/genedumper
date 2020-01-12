#!/usr/bin/env python
#script designed to use blast.db and get accession numbers to pull down.
#This pulls down those with only one choice and pulls longest for those with multiple

def resolve_seqs(blastdb):
    import os, sys, sqlite3, re, time, itertools
    from Bio import Seq, Entrez, SeqIO
    from blastlib.clean_seq_funcs import resolve_seqs, alignment_reg, alignment_comp, alignment_rev_comp, identity_calc, tiling
    from cleanlib.databasing import get_seqs_from_sqldb_GI, get_seqs_from_sqldb_GI_no_gene
    conn = sqlite3.connect(blastdb)
    c = conn.cursor()
    GI_nums_all = set()
    GI_nums_single = set()
    GI_nums_single_GI = []
    mito = set()
    genes = set()
    dic = {}
    dic_ind = {}
    dic_single = {}
    dic_mult = {}
    count2 = 1
    GI_mito_GI = []
    for iter in c.execute("SELECT Gene_name from blast GROUP BY Gene_name;"):
        genes.add(iter[0]) 
    #this gets list of all taxa/genes regardless if they have multiple gene choices or not
    for iter in c.execute("SELECT tc_id, Gene_name, GI, hit_length FROM blast WHERE tc_id != '0' GROUP BY tc_id, Gene_name, GI;"):
        GI_nums_all.add(str(iter[0])+"_"+str(iter[1])+"|"+str(iter[2])+"_"+str(iter[3]))
    #this gets a list of all taxa/genes if they only have one gene choice
    for iter in c.execute("SELECT tc_id, Gene_name, GI, hit_length FROM blast WHERE tc_id != '0' GROUP BY tc_id, Gene_name HAVING COUNT(*) =1;"):
        GI_nums_single.add(str(iter[0])+"_"+str(iter[1])+"|"+str(iter[2])+"_"+str(iter[3]))
    #this give me all the tc_ids that have multiple gene choices
    #tc_id_gene|GI_hit_length
    GI_nums = GI_nums_all-GI_nums_single
    GI_nums_single_GIs = []
    
    #deal with singletones
    for i in GI_nums_single:
        if int(i.split("_")[-1]) > 5000:
            mito.add(i)
    #have to pull out mito/chloro ones separately
    for i in mito:
        GI_nums_single.remove(i)
        GI_mito_GI.append(re.split('_|\|', i)[-2])   
    GI_mito_GI_str = str(GI_mito_GI).replace("[", "(").replace("]", ")")    
    c.execute("UPDATE blast SET Decision='Mito or chloro sequence/Chosen' WHERE GI IN " + GI_mito_GI_str + ";")
    #write singletons
    for i in GI_nums_single:
        GI_nums_single_GIs.append(re.split('_|\|', i)[-2])
    GI_nums_single_GIs_str = str(GI_nums_single_GIs).replace("[", "(").replace("]", ")")
    c.execute("UPDATE blast SET Decision='Only choice/chosen' WHERE GI IN " + GI_nums_single_GIs_str + ";")
    #makes a dictionary of lists for each multiple taxa/gene choice
    #num_gene: ['gi_len']
    for i in GI_nums:
        if i.split("|")[0] in dic.keys():
            dic_list = dic[i.split("|")[0]]
            dic_list.append(i.split("|")[1])
            dic[i.split("|")[0]] = dic_list
        else:
            dic[i.split("|")[0]] = [i.split("|")[1]]    
    
    countall = 0
    #deal with multiples and try and resolve
    print("Trying to resolve sequences")
    for i in dic:
        countall += 1
        #print(i)
        #print(countall)
        print(str(round(float(countall)/float(len(dic))*100, 2))+'%')
        lengths = [int(m.split('_')[1]) for m in dic[i]]
        individual = [dic[i][x] for x, l in enumerate(lengths) if l < 5000]
        mito = [dic[i][x] for x, l in enumerate(lengths) if l >= 5000]
        if len(mito) > 0:
            #will pick the first one if there are multiple ones
            mitoinGI = [mito[0].split("_")[0]]
            c.execute("UPDATE blast SET Decision='Mito or chloro sequence/Chosen' WHERE GI='" + mitoinGI[0] + "';")
            GIS_not_picked_mito = list(set([x.split("_")[0] for x in mito])-set([mitoinGI[0]]))
            if len(GIS_not_picked_mito) != 0:
                GIS_not_picked_mito_str = str(GIS_not_picked_mito).replace("[", "(").replace("]", ")")
                c.execute("UPDATE blast SET Decision='Mito or chloro sequence/Randomly not chosen' WHERE GI IN " + GIS_not_picked_mito_str + ";")
                
            GIS_not_picked_mito = list(set([x.split('_')[0] for x in individual])-set([mitoinGI[0]]))
            if len(GIS_not_picked_mito) != 0:
                GIS_not_picked_mito_str = str(GIS_not_picked_mito).replace("[", "(").replace("]", ")")
                c.execute("UPDATE blast SET Decision='Short or less info/Not chosen' WHERE GI IN " + GIS_not_picked_mito_str + ";")

        elif len(individual) > 0:
            #this should never happen
            # if len(individual) == 1:
            #     c.execute("UPDATE blast SET Decision='Not Chosen' WHERE Name_num='" + individual[0].split("_")[0] + "';")
                
            # else:
#                print(individual)
            result = []
            ranges = {}
            current_start = -1
            current_stop = -1
            whole_length = 0
            #do tiling
            gene = i.split("_", 1)[1]
            for m in [x.split('_')[0] for x in individual]:
                idens, start_stop = tiling([m], gene, blastdb, c)
                start, end = start_stop[0] 
                for iden in idens:
                    if iden > 70:
                        #make dic of lists
                        if (start, end) in ranges.keys():
                            ranges_dic_list = ranges[(start, end)]
                            ranges_dic_list.append(m)
                            ranges[(start, end)] = ranges_dic_list
                        else:
                            ranges[(start, end)] = [m]
                    else:
                        print('Alignment below 70')
                        c.execute("UPDATE blast SET Decision='Sequence does not align well (<70%) to input query sequence/not chosen' WHERE GI='" + m + "';")

                        #print(GIs_to_align)
            if len(ranges) == 0:
                c.execute("UPDATE blast SET Decision='Sequence does not align well (<70%) to input query sequence/not chosen' WHERE Name_num='" + i + "';")
                print('All alignments below 70, printing to file')
                with open("fulllength_hand_check.txt", "a") as a:
                    a.write(i + '\n')
            else:
                #get merged range
                for start, stop in sorted(ranges.keys()):
                    if start > current_stop:
                        result.append((start, stop))
                        current_start, current_stop = start, stop
                    else:
                        current_stop = max(current_stop, stop)
                        result[-1] = (current_start, current_stop)
                for n in result:
                    whole_length += n[1] - n[0] + 1
                #go through each combination of ranges to get 95% of whole range
                for L in range(1, len(ranges)+1):
                    max_perc = 0
                    for subset in itertools.combinations(ranges.keys(), L):
                        comb_length = 0
                        #for each combination, get merged length
                        current_start = -1
                        current_stop = -1
                        result = []
                        for start, stop in sorted(subset):
                            if start > current_stop:
                                result.append((start, stop))
                                current_start, current_stop = start, stop
                            else:
                                current_stop = max(current_stop, stop)
                                result[-1] = (current_start, current_stop)
                        for x in result:
                            comb_length += x[1] - x[0] + 1
                        if whole_length >= comb_length:
                            perc = float(comb_length)/float(whole_length)
                            if perc > max_perc:
                                max_perc = perc
                                max_subset = set()
                                max_subset.add(subset)
                            elif perc == max_perc:
                                max_subset.add(subset)
                            else:
                                pass  
                        else:
                            pass
                    # goes through all combinations in a level before breaking
                    if max_perc > .95:
                        break
                final_tiling = [(0, 0)]*L
                for combination in max_subset:
                    for x, comb_frag in enumerate(sorted(combination)):
                        if comb_frag[1] - comb_frag[0] + 1 > final_tiling[x][1] - final_tiling[x][0] + 1:
                            final_tiling[x] = comb_frag
                if (0, 0) in final_tiling:
                    possible_GIs = []
                else:
                    possible_GIs = [ranges[x] for x in final_tiling]
                GIS_not_picked_tiling = list(set([x.split('_')[0] for x in individual])-set([GI for GI_tiles in possible_GIs for GI in GI_tiles]))
                if len(GIS_not_picked_tiling) != 0:
                    GIS_not_picked_tiling_str = str(GIS_not_picked_tiling).replace("[", "(").replace("]", ")")
                    c.execute("UPDATE blast SET Decision='Better tiling/Not chosen' WHERE GI IN " + GIS_not_picked_tiling_str + ";")

                count = 1
                chosen_GIs = []
                for m in possible_GIs:
                    if len(m) == 1:
                        c.execute("UPDATE blast SET Decision='Only or best choice in tiling analysis/chosen' WHERE GI='" + m[0] + "';")
                    else:
                        GIs_to_pick = resolve_seqs(m, blastdb, gene, c)
                        GIs_to_pick_str = str(GIs_to_pick).replace("[", "(").replace("]", ")")
                        if len(GIs_to_pick) == 1:
                            c.execute("UPDATE blast SET Decision='Longest or most info (not checked), tile "+str(count)+" /Chosen' WHERE GI IN " + GIs_to_pick_str + ";")
                            GIS_not_picked = list(set(m)-set(GIs_to_pick))
                            chosen_GIs.append(GIs_to_pick)
                            GI_not_picked_str = str(GIS_not_picked).replace("[", "(").replace("]", ")")
                            c.execute("UPDATE blast SET Decision='Short or less info, tile "+str(count)+"/Not chosen' WHERE GI IN " + GI_not_picked_str + ";")
                            count += 1
                        else:
                            c.execute("UPDATE blast SET Decision='To cluster analysis/Chosen' WHERE GI IN " + GIs_to_pick_str + ";")
                            GIS_not_picked = list(set(m)-set(GIs_to_pick))
                            chosen_GIs.append(GIs_to_pick)
                            if len(GIS_not_picked) != 0:
                                GI_not_picked_str = str(GIS_not_picked).replace("[", "(").replace("]", ")")
                                c.execute("UPDATE blast SET Decision='Short or less info/Not chosen' WHERE GI IN " + GI_not_picked_str + ";")
                        #merge ranges - if same number of ranges as original, keep all,
                        #want the least number to overlap 95% of whole range
                        #try each comb from low to high and if hits 95%, choose those
                        #if multiple higher than 95%, choose one with best %
                        #if multiple with same % and same numb of combs- multiple
                
    #print(dic)

    conn.commit()
    conn.close()
