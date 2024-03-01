This is a modified version of the POEM pipeline with limited functionality, which can only predict operons in genome assemblies. This version was created for predicting operons in Comparative Genomic Content Management System. It can neither train operon models nor work with metagenomes and short reads. Full version is available at [https://github.com/Rinoahu/POEM_py3k]


## Requirements


This pipeline is available on Linux systems. Make sure that you have the following installed on linux

1.  [Anaconda](https://www.anaconda.com/ "https://www.anaconda.com/") for Python 3.7
2.  [conda](https://conda.io/en/latest/ "https://conda.io/en/latest/")

    make sure to add path of conda to $PATH environment variable

## Linux Installation

Installation is simple once Anaconda and Conda are installed. Type or paste the following commands into your terminal in whichever subfolder you want to keep POEM.

```
$ git clone https://github.com/aekazakov/POEM_py3k

$ cd ./POEM_py3k

$ bash ./install.sh <condaenv>
```
where condaenv is a desired name for conda environment.

The installation script calls conda to install all the necessary python packages and software, as well as the COG database. To avoid conflict, it will create a new conda environment.

## Testing the Installation

The example directory contains genome fasta and protein fasta files, run  ```runme.sh``` to test the pipeline.

To check if POEM is installed correctly, the users can make a quick and simple test:
```
$ cd ./example

$ bash ./runme.sh
```

This should output a file called input.fsa.adjacency containing 41 lines.


## Usage
This POEM version is recommended for finding operons in assembled genomes only. It requires two files: genome fasta file and protein fasta file.

Run the bin/run_poem_cgcms.sh script.

## Output


POEM will create several files in the directory with input files, but you only need one of them:

        input.fsa.adjacency:
		This file is the predicted operonic adjacent genes. This file contains 11 columns: col 1 is the gene 1's id; col 2 is the chromosome|scaffold id where gene 1 is located; col 3-5 are strand, start and end of gene 1; col 6 is gene 2's id; col 7 is the chromsome|scaffold id where gene 2 is located; col 8-11 are strand, start and end of gene 2. For example:

		CPIBHFJP_00001|Prodigal:2.6|772_aa|+|255|857$$gi|985000614|gb|CP014225.1|       gi|985000614|gb|CP014225.1|     +       255     857     CPIBHFJP_00002|Prodigal:2.6|1014_aa|+|883|1308$$gi|985000614|gb|CP014225.1|     gi|985000614|gb|CP014225.1|     +       883     1308    False
		CPIBHFJP_00002|Prodigal:2.6|1014_aa|+|883|1308$$gi|985000614|gb|CP014225.1|     gi|985000614|gb|CP014225.1|     +       883     1308    CPIBHFJP_00003|Prodigal:2.6|1165_aa|+|1586|1693$$gi|985000614|gb|CP014225.1|    gi|985000614|gb|CP014225.1|     +       1586    1693    False
		...


