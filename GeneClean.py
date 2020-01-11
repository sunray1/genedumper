#!/usr/bin/env python

import sys, os
from cleanlib.multiple import resolve_seqs
from cleanlib.blast_sp_parse import test_resolved_seqs, first_blast
from cleanlib.getseqsfromGenbank import pullseqs
from cleanlib.cluster_analysis import cluster
from argparse import ArgumentParser, RawTextHelpFormatter

argp = ArgumentParser(description='Takes sqlite database with sequences resolved by tc_id and tries to choose the best sequence for each tcid/gene pair\nmodule load python/2.7.6 muscle', formatter_class=RawTextHelpFormatter)
argp.add_argument('-s', '--steps', help='the steps the program will run \n0 - Run initial resolver\n1 - Pull sequences for BLAST check\n2 - Check resolved sequences\n3 - Cluster Analysis\n4 - Pull cleaned sequences down\n')
argp.add_argument('-b', '--blastdb', help='the name of the sqlite database (blast_results.db by default)')
argp.add_argument('-t', '--taxdb', help='the name of the taxonomy database', required=True)
argp.add_argument('-e', '--email', help='user email used for NCBI')
argp.set_defaults(blastdb='blast_results.db', taxdb=False)
args = argp.parse_args()
blastdb = args.blastdb

if len(sys.argv) == 1:
    argp.print_help()
    sys.exit()

if not args.steps:
    steps = '01234'
else:
    steps = args.steps
if not os.path.isfile(args.blastdb):
    sys.exit("Error: no blast sqlite file called " + blastdb)
if not args.taxdb:
    sys.exit('Error: need a taxonomy database')
elif args.taxdb:
    if not os.path.isfile(args.taxdb):
        sys.exit("Error: no taxonomy file called " + args.taxdb)
if '4' in steps and not args.email:
    sys.exit("Error: NCBI needs email -e to pull down sequences")
 
 
taxdb = args.taxdb
#run multiple.py - breaks everything apart and tries to initially resolve it


if '0' in steps:
    print('Running initial resolver')
    resolve_seqs(blastdb, email)

if '1' in steps:
    first_blast(blastdb)

if '2' in steps:
    list_in = [f for f in os.listdir(".") if f.endswith("_qseqs.fa.xml")]
    # print(list_in)
    print('Checking resolved sequences')
    for f in list_in:
        test_resolved_seqs(f, blastdb, taxdb, email)
        
if '3' in steps:
    print('Running cluster analysis for unresolved choices')
    cluster(blastdb, taxdb, email)
    
if '4' in steps:
    email = args.email
    print("Pulling down final cleaned sequences")
    pullseqs(blastdb, email)
