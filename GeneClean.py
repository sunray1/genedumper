#!/usr/bin/env python

import sys, os
from cleanlib.multiple import resolve_seqs
from cleanlib.blast_sp_parse import test_resolved_seqs
from cleanlib.getseqsfromGenbank import pullseqs
from cleanlib.cluster_analysis import cluster
from argparse import ArgumentParser, RawTextHelpFormatter

argp = ArgumentParser(description='Takes sqlite database with sequences resolved by tc_id and tries to choose the best sequence for each tcid/gene pair\nmodule load python/2.7.6 muscle', formatter_class=RawTextHelpFormatter)
argp.add_argument('-b', '--blastdb', help='the name of the sqlite database (blast_results.db by default)')
argp.add_argument('-t', '--taxdb', help='the name of the taxonomy database')
argp.set_defaults(blastdb='blast_results.db', taxdb=False)
args = argp.parse_args()
blastdb = args.blastdb

if not os.path.isfile(args.blastdb):
    sys.exit("Error: no blast sqlite file called " + blastdb)
if not args.taxdb:
    sys.exit('Error: need a taxonomy database')
elif args.taxdb:
    if not os.path.isfile(args.taxdb):
        sys.exit("Error: no taxonomy file called " + args.taxdb)
    
 
taxdb = args.taxdb
#run multiple.py - breaks everything apart and tries to initially resolve it
print('Running initial resolver')
#resolve_seqs(blastdb)

list_in = [f for f in os.listdir(".") if f.endswith("accession_nums_resolved.txt")]

print("Pulling down seqs for blasting")
for f in list_in:
    print("Pulling down " + f)
    pullseqs(f)
#list_in = ["CAT_accession_nums_resolved.fa"]
list_in = [f for f in os.listdir(".") if f.endswith("accession_nums_resolved.fa")]
# print(list_in)
print('Checking resolved sequences')
for f in list_in:
	test_resolved_seqs(f, blastdb, taxdb)

if os.path.getsize("multiple_gene_choices.txt") > 0:
    print('Running cluster analysis for unresolved choices')
    cluster(blastdb, taxdb)
else:
    os.remove("multiple_gene_choices.txt")
