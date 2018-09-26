#!/usr/bin/env python

from Bio import SeqIO
import sys, sqlite3

def number(filein):
    count = 1
    recs_changed = {}
    record_to_print = set()
    records = list(SeqIO.parse(filein, "fasta"))
    for record in records:
        print(float(count)/float(len(records)))
        count+=1
        record.description = ''
        id = record.id
        if id in recs_changed.keys():
            num = int(recs_changed[id].split('_')[-1])
            recs_changed[id] = id+'_'+str(num+1)
            record.id = id+'_'+str(num+1)
            
        else:
            recs_changed[id] = id+'_1'
        record_to_print.add(record)
    outfile = open(filein.split(".")[0] + "_renamed_with_nums.fa", "w")
    SeqIO.write(record_to_print, outfile, "fasta")