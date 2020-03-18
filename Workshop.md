# GeneDumper Workshop
---
# Setup


1. Login to cluster
We are using the orange drive today for space reasons.

```sh
cd /orange/kawahara/GDWorkshop
mkdir $USER
cd $USER
```
2. Load modules

```sh
module load git python3 muscle R
```

3. Download GeneDumper

```sh
git init
git clone https://github.com/sunray1/GeneDumper.git
cd GeneDumper
```
---
# Inputs



Two inputs are needed: 

![Flowchart of GeneDumper](https://github.com/sunray1/Images/blob/master/GeneDumper.jpg?raw=true)

1. Taxonomy file in .csv format for your Taxa of Interest (TOI) (example [here](https://github.com/sunray1/GeneDumper/blob/master/example_files/Butterfly/Lama_Pieridae.csv))
      -  This can be built on your own (ex. Lamas list) or downloaded from various resources (ex. WorldFloraOnline, ITIS)
      -  Alternatively we can download the taxonomy from GenBank
      -  We will convert this to an .sql file

2. Fasta file (example [here]( https://github.com/sunray1/GeneDumper/blob/master/example_files/Butterfly/D_plexippus_probes.fas))
      -  Sequence names can't have spaces 
      -  I usually just name this with the gene name (ie. COI, Ef1a)


### Taxonomy file

You should have your taxa of interest (preferably a family or something smaller to make this workshop go quickly) that is a part of [Genbank's taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy).

##### **Getting an API Key**

In order to do this, each of you will also need a GenBank API Key. This allows you to ping GenBank's API with a much higher frequency. Luckily, you can log in using your GatorLink account (if you don't already have one). 

https://www.ncbi.nlm.nih.gov/account/register/

![NCBI GatorLink Registration](https://github.com/sunray1/Images/blob/master/ncbi.JPG?raw=TRUE)

Once you have logged in, click on your name in the top right corner to go to your account settings and scroll down to API Key Management. Click on `Create an API Key` to generate your own.

![API Key Management](https://github.com/sunray1/Images/blob/master/api.JPG?raw=TRUE)

##### **Downloading your taxonomy using R**
 
 I've made a taxonomy folder availiable inside the GeneDumper folder for people to work out of called `taxonomy_files`. This file contains scripts to create your own taxonomy (like we're doing here) and libraries to convert it into an SQL database.
```sh
$ cd taxonomy_files
$ ls
empty.db      load_taxonomy.py  taxo_of_taxa_ex.R
example.conf  taxolib           taxonomy.R    taxo_of_taxa.R

```

1. Open `taxo_of_taxa.R` in a text editor (I prefer to use vim) so we can edit some parameters.
    - Uncomment line so you can download the updated version of taxize
    - Add your API Key in
    - Add your Taxon of Interest in
    - You can add or remove Ranks of Interest here (don't put Genus or Species, these are automatic) (Don't put subgenus for now, theres an error I need to fix)
    - Add your division filter

![taxo_of_taxa.R edits](https://github.com/sunray1/Images/blob/master/taxoR.PNG?raw=TRUE)

The division filter is useful when there are multiple rankings with the same name. This doesn't really matter if we were using RStudio, but here we want to add it to reduce interactiveness. You can find this on [Genbank's taxonomy site](https://www.ncbi.nlm.nih.gov/taxonomy). If you search for your taxon of interest, it is the NCBI BLAST name, **NOT** the Genbank common name.

[Example](https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=9787&lvl=3&lin=f&keep=1&srchmode=1&unlock):
![Division filter](https://github.com/sunray1/Images/blob/master/divison_filter.PNG?raw=TRUE)


2. Once your file is saved, go ahead and run it!

```sh
Rscript taxo_of_taxa.R
```

Once it has finished, you should have a `Taxonomy.csv` file in your directory.

##### **Converting your Taxonomy into an SQL file**

Inside this folder is an empty SQL database called `empty.db` that we need to populate with our data. To do so, we need to edit the `example.config` file. This is where we tell the database what our ranks are and if we have synonomys. 

1. Open `example.config` in your text editor so we can edit some parameters

Working from the top down, there are six sections and a few lines we need to change. 

| Section Title        | Line to change | Notes|  
| :-------------: |:-------------:|:--------|
| [main]      | inputcsv = | Taxonomy.csv
| [main]      | ranksys = | Animals is 1 & if you get plant data off GenBank, it will be 2
| [taxonomy]      |  name =      | Can be whatever
| [root] | rank =       | The rank of your TOI
| [root] | name = | The name of your TOI|
| [synonyms] | colname = Synonyms | Comment out
| [synonyms] | separator = ,| Comment out
|[rank_mappings] | | This should be the ranks you picked in the R script

2. Once the config file is saved, run the code below to populate the SQL file.
```sh
python load_taxonomy.py -l none -d empty.db example.conf
```

3. And congrats, now you have a SQL database containing your data. Go ahead and rename it if you'd like (since it's not empty anymore)

```sh
mv empty.db Perissodactyla.db
```

### FASTA file

In order to run, GeneDumper needs a FASTA file containing one sequence to use as a seed sequence when searching. Ideally you'd use full-length, high quality sequences here. You can determine which sequence to use by looking at the literature to see what people use or by looking at the sequences themselves.

1. Go to your TOI's [taxon site](https://www.ncbi.nlm.nih.gov/taxonomy) and click on the Nucleotide link. 

![Nucleotide record link](https://github.com/sunray1/Images/blob/master/ntrecs.PNG?raw=TRUE)

I decided to go with BRCA1 and 12S due to these two papers:
https://link.springer.com/article/10.1007/s002399910002
https://onlinelibrary.wiley.com/doi/10.1111/j.1096-3642.2011.00752.x

2. Grab a sequence from a respresentative of your TOI and paste it into a new file using a text editor. I like to rename these sequences so they make more sense - be sure to remove spaces.

---
# GeneDump
Now we can actually start running GeneDumper!

GeneDump is the first part that we are going to run. This is the part that does the initial blast (takes the longest) and name resolution of your hits to your taxonomy.

```sh
$ python GeneDump.py -h
usage: GeneDump.py [-h] [-b BLASTDB] [-s STEPS] [-t TAXDB] [-c CALCSTATS]
                   [-f FASTAIN] [-l] [-e EMAIL]

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
  -l, --html            calculates an output html file with hero stats ab
  -e EMAIL, --email EMAIL
                        email for NCBI
```

There are seven steps you will have to sequentially run. Lets just run the blast step for now (0). The amount of time this part takes can vary widely depending on how many sequences you're searching, how many people are using NCBI and how large your TOI is.

```sh
python3 GeneDump.py -s 0 -t taxonomy_files/Perissodactyla.db -f seedseqs.fa
```

Now lets go ahead and run the rest of the steps. These shouldn't take too long.

```sh
python3 GeneDump.py -s 123456 -t taxonomy_files/Perissodactyla.db -e sunray1@ufl.edu
```

### Output files
There are a few output files here that you can look at to see what changed in our database.

1. `blast_results.db` - this is the blast database containing all returned sequences and their metadata. We use this as input into the next step.
2. `spell.txt`, `ncbi.txt`, `epithet.txt` - these are files you can look at to see if there were any misspellings or incorrectly labeled sequences. There shouldn't be much in these (maybe some in `ncbi.txt`) since we used a GenBank taxonomy.

### Fun things
There are some fun stats we can calculate at this point if we want to. 

 - To get a presense/absense matrix of our data, we can use the -c command. This can be used at any taxonomic level, but perhaps species is the most important.

```sh
python3 GeneDump.py -t taxonomy_files/Perissodactyla.db -c Species
```

 - To generate an html file containing some hero stats, we can use the -l command. 

```sh
python3 GeneDump.py -t taxonomy_files/Perissodactyla.db -l
```

---
# GeneClean

Cleaning the sequences is also fairly easy to do. 

```sh
$ python GeneClean.py -h
usage: GeneClean.py [-h] [-s STEPS] [-b BLASTDB] -t TAXDB [-e EMAIL]

Takes sqlite database with sequences resolved by tc_id and tries to choose the best sequence for each tcid/gene pair
module load python/2.7.6 muscle

optional arguments:
  -h, --help            show this help message and exit
  -s STEPS, --steps STEPS
                        the steps the program will run
                        0 - Run initial resolver
                        1 - Pull sequences for BLAST check
                        2 - Check resolved sequences
                        3 - Cluster Analysis
                        4 - Pull cleaned sequences down
  -b BLASTDB, --blastdb BLASTDB
                        the name of the sqlite database (blast_results.db by default)
  -t TAXDB, --taxdb TAXDB
                        the name of the taxonomy database
  -e EMAIL, --email EMAIL
                        user email used for NCBI
```

It is structured the same way as GeneDump, in steps. This step choses 'best' sequences based on the decision flowchart below.

image

Lets run this.

```sh
python GeneClean.py -t taxonomy_files/Perissodactyla.db -e sunray1@ufl.edu
```

We can also re-run the html command again to get an idea of our decisions.

And we are finished! We have cleaned sequences :)

---

I made a quick, fun little tree based on my sequences for kicks and giggles - 

![FastTree phylogeny](https://github.com/sunray1/Images/blob/master/tree.png?raw=TRUE)
