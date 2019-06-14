# GeneDumper

GeneDumper is an auto-updating software tool to clean and maintain GenBank data for a user’s taxa and genes of interest. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. GeneDumper is written in [Python](https://www.python.org/)(2.7+), uses [BioPython](https://biopython.org/) and [SQLite](https://www.sqlite.org/index.html) libraries and command line interface, and runs on a GNU/Linux operating system. We welcome feedback, bug reports, or code contributions from users.

### How it Works

Required inputs:
* Taxonomy in .csv format (example)
1. Taxonomy Databasing - 
2. Initial BLAST and Species Name Resolution
3. Cleaning and Validation

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
git init
git clone https://github.com/sunray1/GeneDumper.git
```


## Usage and Arguments



### Initial BLAST and SQL database formation

Explain what these tests test and why

```
Give an example
```

### Sequence Cleaning and Validation

Explain what these tests test and why

```
Give an example
```



## Versioning
#### July 2019: Version 0.8.0
Initial Release

## Authors

* **Chandra Earl** - *Primary Author* - [Github](https://github.com/sunray1)
* **Brian Stucky** - *Taxonomy databasing* - [Github](https://github.com/stuckyb)
* **Rob Guralnick** - *Principal Investigator* - [Github](https://github.com/robgur)
* **Akito Kawahara** - *Principal Investigator* - [Github]()

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
