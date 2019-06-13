
# First Python Notebook

Hello World


```python
import os, sys
```


```python
ls
```

    [0m[34;42malignlib[0m/  [34;42mexample_files[0m/  [01;32mREADME.md[0m*       [01;32muntitled1.txt[0m*
    [34;42mblastlib[0m/  [01;32m__init__.py[0m*    [01;32mrun_scripts.py[0m*  [01;32muntitled.txt[0m*
    [34;42mcleanlib[0m/  [01;32mReadme.ipynb[0m*   [34;42mtaxonomy_files[0m/



```python
%run -i run_scripts.py -h
```

    usage: run_scripts.py [-h] [-b BLASTDB] [-s STEPS] [-t TAXDB] [-c CALCSTATS]
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

