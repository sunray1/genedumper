# GeneDumper

GeneDumper is an auto-updating software tool to clean and maintain GenBank data for a user’s taxa and genes of interest. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. GeneDumper is written in [Python](https://www.python.org/)(2.7+), uses [BioPython](https://biopython.org/) and [SQLite](https://www.sqlite.org/index.html) libraries and command line interface, and runs on a GNU/Linux operating system. We welcome feedback, bug reports, or code contributions from users.

### How it Works

Required inputs for each step:
1. Taxonomy Databasing
    a. Taxonomy in .csv format (can contain synonyms in separate column) OR
    b. List of Species (can contain synonyms in separate column)
    c. .config file describing taxonomy (can edit given one)
2. Initial BLAST and Species Name Resolution
    a. One FASTA file containing one of each of the loci of interest. Ideally these sequences need to be full length.
3. Cleaning and Validation
    a. SQL database of BLAST hits created by step 2.

Example inputs:
    [Taxonomy in csv format](https://github.com/sunray1/GeneDumper/blob/master/example_files/Butterfly/Butterflynet2tax_syns.csv)
    [FASTA file](https://github.com/sunray1/GeneDumper/blob/master/example_files/Butterfly/D_plexippus_probes.fas)
    [.config file](https://github.com/sunray1/GeneDumper/blob/master/example_files/Butterfly/ButterflyNet.config)


### Prerequisites

[SQLite](https://www.sqlite.org/download.html) is used for structuring and querying of the database itself, using Python as a wrapper for the user interface. BioPython is used for BLAST and calculating alignments. It is designed to run in a command-line environment and assumes all prerequisites are in the user's PATH. While the both the sqlite3 command line tool and SQLite python libraries needed, the libraries are built into python2.5+.


To check that your prerequisites are installed and running correctly, make sure the following commands finish without error:
```
$ sqlite3
SQLite version 3.11.0 2016-02-15 17:29:24
Enter ".help" for usage hints.
sqlite>
```
```
$python
Python 2.7.12 (default, Nov 19 2016, 06:48:10)
[GCC 5.4.0 20160609] on linux2
Type "help", "copyright", "credits" or "license" for more information.
>>> import sqlite3
>>> import Bio
```


### Installing

Currently, scripts only need to be pulled directly from Github and given executable permissions: ```chmod +x GeneDumper/```. Packages will be availiable soon for pip installaton.

Using command line git-clone:
```
$ git init
$ git clone https://github.com/sunray1/GeneDumper.git
```


## Usage and Arguments

### Taxonomy database formation

The goal of the first step is to get your taxonomy for your taxa of interest into a SQL database. This can be done using a given taxonomy (in .csv format) or the software can infer the taxonomy for you based on NCBI's taxonomy. 

If you need an inferred taxonomy given a list of species/synonyms:

```
Give an example
```

Once you have a taxonomy OR are providing one you need two other input files. Both of these files are included in the taxonomy_files folder.
1. .config file (can edit the one that is included)
2. empty taxonomy databse (included)

The .config file contains all options needed by load_taxonomy.py to parse the inputted .csv taxonomy file. These options will change depending on the layout of the .csv - see example_files for further help.

The empty taxonomy database (empty.db) is a structured SQL database, without the data.
 
The script to load your taxonomy into the empty database is also located in the taxonomy_files folder and runs as follows (database and .conf file names can be changed according to your project):

```
$ python load_taxonomy.py -l none -d empty.db example.conf
```
It uses libraries found in the taxolib folder and should be run in the same location as this folder. This will load your data from your .csv file (location noted in .conf file) into empty.db. This database will be used in subsequent steps as the main taxonomy.

### Initial BLAST and SQL database formation

The script to run the intial BLAST and put results into a SLQ database is called run_scripts.py. It uses libraries in the blastlib folder and should be run in the same location as this folder. It runs certain steps based on user input. -b is the name of the outputted BLAST sequence SQL database, while -t is the inputted taxonomy database made in the step above. 

A FASTA file containing your loci of interest is required for the BLAST step.

```
$ python run_scripts.py -h
usage: GeneDump.py [-h] [-b BLASTDB] [-s STEPS] [-t TAXDB] [-c CALCSTATS]
                      [-f FASTAIN]

Takes .xml files from blast, filters seqs and resolves names

optional arguments:
  -h, --help            show this help message and exit
  -b BLASTDB, --blastdb BLASTDB
                        the name of the sqlite database (blast_results.db by default)
  -s STEPS, --steps STEPS
                        the steps the program will run
                        0 - run blast using -f fasta file
                        1 - .xml to .csv
                        2 - .csv to .sql
                        3 - resolve names
                        4 - check ncbi for correct names
                        5 - check epithets for matches
                        6 - check spelling errors
                        default is to run all steps
  -t TAXDB, --taxdb TAXDB
                        the name of the taxonomy database
  -c CALCSTATS, --calcstats CALCSTATS
                        calculate statistics based on the taxonomy level given, ie Species, Genus
  -f FASTAIN, --fastain FASTAIN
                        the fasta file that blast will use

```

An example to run all the steps in this script:

```
$ python GeneDump.py -b output_blast.db -t input_taxonomy.db -f FASTA_file -s 0123456
```
Alternatively to run all steps, the -s arguement can be ignored as the script runs all steps by default.

The -c option allows the user to produce a presence/absense matrix depicting what was found in the initial blast. This is usually at the Species level, but other levels can be specified instead (Genus, Family, etc.).



### Sequence Cleaning and Validation

This last step is useful if you wish to pull down a cleaned version of your sequences. The script for this is GeneClean.py and uses the alignlib and cleanlib libraries and should be run in the same location as these folders. This step can take ~24 hours depending on your blast database size.

This step only takes two inputs - the blast database you wish to pull from and the taxonomy database associated with the data. 

```
$ python GeneClean.py -b blast_db -t taxonomy_db
```

## Outputs

ouputs and files made from each step

## Versioning
#### July 2019: Version 0.8.0
Initial Release

## Authors

* **Chandra Earl** - *Primary Author* - [Github](https://github.com/sunray1)
* **Brian Stucky** - *Taxonomy databasing* - [Github](https://github.com/stuckyb)
* **Rob Guralnick** - *Principal Investigator* - [Github](https://github.com/robgur)
* **Akito Kawahara** - *Principal Investigator* - [Github]()

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Step in CSV to SQL step of GeneDump forked from [Rufus Pollock](https://github.com/rufuspollock/csv2sqlite). Many thanks!
