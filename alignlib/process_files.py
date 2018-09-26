#!/usr/bin/env python

import os, sys, sqlite3
#from cleanlib.getseqsfromGenbank import pullseqs
def process(files, blastdb, pullseqs):
    conn = sqlite3.connect(blastdb)
    fasta_files = []
    GI_files = []
    c = conn.cursor()
    count = 1
    dic = {}
    
    for f in files:
        with open(f) as o:
            line = o.readline()
            if line[0] == ">":
                fasta_files.append(f)
            else:
                GI_files.append(f)
    cat_command = "cat " + "%s "*len(GI_files) + "> allGIs.txt"
    os.system(cat_command % tuple(GI_files))
    
    with open("allGIs.txt") as o:
        line = o.readline()
        while line:
            print(float(count)/float(6))
            for iter in c.execute("SELECT Gene_name FROM blast WHERE GI='"+ line.strip() + "'"):
                gene = str(iter[0])
            if gene in dic.keys():
                dic_list = dic[gene]
                dic_list.append(line.strip())
                dic[gene] = dic_list
            else:
                dic[gene] = [line.strip()]
            line = o.readline()
            count += 1
    filelist = []
    for g in dic:
        filelist.append(g+"_GIs.txt")
        with open(g+"_GIs.txt", "w") as o:
            for m in dic[g]:
                o.write(m+'\n')
    
    for f in filelist:
        pullseqs(f)
    
    
    cat_command = "cat " + "%s "*len(fasta_files) + "COI_trnL_COII_GIs.fa > COI_trnL_COII_GIs_all.fa"
    os.system(cat_command % tuple(fasta_files))
    os.remove("allGIs.txt")