#!/usr/bin/env python

#Script designed as a 'command center' to parse arguments, check for basic errors and run all other functions. One function it does include
#is the blast function however as this is a one line command

import sys, os
from blastlib.edit_csv import edit
from blastlib.csv2sqlite import convert
from blastlib.name_res import find_id
from blastlib.check_ncbi import ncbi
from blastlib.epithet_check import epithet
from blastlib.us_to_a import spelling
from blastlib.stats import get_stats
from argparse import ArgumentParser, RawTextHelpFormatter

argp = ArgumentParser(description='Takes .xml files from blast, filters seqs and resolves names', formatter_class=RawTextHelpFormatter)
argp.add_argument('-b', '--blastdb', help='the name of the sqlite database (blast_results.db by default)')
argp.add_argument('-s', '--steps', help='the steps the program will run \n0 - run blast using -f fasta file\n1 - .xml to .csv \n2 - .csv to .sql \n3 - resolve names \n4 - check ncbi for correct names \n5 - check epithets for matches \n6 - check spelling errors\ndefault is to run all steps\n')
argp.add_argument('-t', '--taxdb', help='the name of the taxonomy database')
argp.add_argument('-c', '--calcstats', help='calculate statistics based on the taxonomy level given, ie Species, Genus')
argp.add_argument('-f', '--fastain', help='the fasta file that blast will use')
argp.set_defaults(blastdb='blast_results.db', taxdb=False)
args = argp.parse_args()
blastdb = args.blastdb
steps = args.steps

if len(sys.argv) == 1:
    argp.print_help()
    sys.exit()

if not args.calcstats and not args.steps:
    steps = '0123456'
elif args.calcstats and not args.steps:
    steps = ''
#Error calling
if not args.taxdb and ('3' in steps or '4' in steps or '5' in steps or '6' in steps or args.calcstats):
    sys.exit('Error: need a taxonomy database')
elif args.taxdb:
    if not os.path.isfile(args.taxdb):
        sys.exit("Error: no taxonomy file called " + args.taxdb)
    taxdb = args.taxdb
if not os.path.isfile(blastdb) and '2' not in steps and ('3' in steps or '4' in steps or '5' in steps or '6' in steps):
    sys.exit("Error: no blast sqlite file called " + blastdb)
    
#Start running functions
if '0' in steps:
    print('Running step 0: run blast\n')
    if not args.fastain:
        sys.exit("Error: needs a fasta file")
    else:
        if not os.path.isfile(args.fastain):
            sys.exit('Error: no fasta file called ' + args.fastain)
        from Bio.Blast import NCBIWWW
        from Bio import SeqIO
        records = SeqIO.parse(args.fastain, format="fasta")
        for record in records:
            id = record.id
            print("Blasting " + id)
            result_handle = NCBIWWW.qblast("blastn", "nt", record.seq, entrez_query='((Papilionoidea[Organism]) OR Hedylidae[Organism]) OR Hesperiidae[Organism]', word_size=7, hitlist_size=50000)
            with open(id+".xml", "w") as save_file:
                save_file.write(result_handle.read())
                result_handle.close()
    print("Step 0 completed\n--------------------")
if '1' in steps:
    print("Running step 1: converting .xml to .csv.....\n")
    xmlfiles = [f for f in os.listdir('.') if f.endswith('.xml')]
    if len(xmlfiles) == 0:
        sys.exit("Error: no .xml files found")
    for file in xmlfiles:
        print("Converting " + file + "...")
        test = edit(file)
        test.xmltocsv()
        print(file + " converted\n")
    print("Step 1 completed\n----------------------")
if '2' in steps:
    print("Running step 2: converting .csv to .sql......\n")
    if os.path.isfile(blastdb):
        os.remove(blastdb)
    csvfiles = [f for f in os.listdir('.') if f.endswith('.csv')]
    if len(csvfiles) == 0:
        sys.exit("Error: no .csv files found")
    for file in csvfiles:
        print("Converting " + file + "...")
        convert(file, blastdb, 'blast')
        print(file + " converted\n")
    print("Step 2 completed\n----------------------")
if '3' in steps:
    print("Running step 3: resolving names......\n")
    find_id(taxdb, blastdb)
    print("Step 3 complete\n----------------------")
if '4' in steps:
    print("Running step 4: checking ncbi for names.....\n")
    ncbi(taxdb, blastdb)
    find_id(taxdb, blastdb)
    print("Step 4 complete, ncbi.txt file contains changes\n----------------------")
if '5' in steps:
    print("Running step 5: checking epithets against species and subspecies.....\n")
    epithet(taxdb, blastdb)
    find_id(taxdb, blastdb)
    print("Step 5 complete, epithet.txt file contains changes\n----------------------")
if '6' in steps:
    print("Running step 6: checking names for mispellings.....\n")
    spelling(taxdb, blastdb)
    find_id(taxdb, blastdb)
    print("Step 6 complete, spelling.txt file contains changes and changes that must be made by hand\n----------------------")
if args.calcstats:
    level = args.calcstats
    print("Calculating statistics for " + level + " level.....\n")
    get_stats(taxdb, blastdb, level)
    print("Calculating statistics completed\n--------------------")
print('Done')