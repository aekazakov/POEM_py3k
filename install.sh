#!/bin/bash
if ! [ $1 ]
then
	echo "Required parameter: conda environment name for installation"
	exit
else
	CONDAENV=$1
fi
CONDAENV=$1
# get the scripts path
SCRIPT=`realpath -s $0`
SCRIPTPATH=`dirname $SCRIPT`

# build a new environment
conda create -y -n $CONDAENV python=3.7

eval "$(conda shell.bash hook)"
conda activate $CONDAENV 
# install third part software and packages:
conda install -y -c bioconda diamond idba cd-hit biopython numpy
conda install -y -c bioconda prokka=1.12 
conda install -y -c bioconda keras=2.2.4 networkx art 
conda install -y -c biobuilds perl=5.22
conda deactivate
# Patch installation
conda activate $CONDAENV 
pip install h5py==2.10.0 --force-reinstall
pip uninstall --yes numpy
pip uninstall --yes scipy
pip uninstall --yes pandas
pip uninstall --yes matplotlib
pip uninstall --yes tensorflow
pip install --no-input numpy==1.16.4
pip install --no-input scipy==1.1
pip install --no-input pandas==0.25.1
pip install --no-input matplotlib==3.4
pip install --no-input tensorflow==1.13.1

# download the cog2014 database.
echo "Download COG database"
wget -q -c ftp://ftp.ncbi.nih.gov/pub/COG/COG2014/data/*2014* -P $SCRIPTPATH/database/cog/cog2014

# unzip cog fasta file and format
echo "Unzip cog fasta and format the file by diamond"
gunzip $SCRIPTPATH/database/cog/cog2014/prot2003-2014.fa.gz
diamond makedb --in $SCRIPTPATH/database/cog/cog2014/prot2003-2014.fa -d $SCRIPTPATH/database/cog/cog2014/prot2003-2014.fa
conda deactivate