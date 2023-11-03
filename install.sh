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

eval "$(conda shell.bash hook)"
# build a new environment
conda create -y -n $CONDAENV python=3.7
conda activate $CONDAENV 
# install third part software and packages:
conda install -y -c bioconda numpy biopython
pip install --no-input tensorflow==1.13.1
pip install --no-input keras==2.2.4
pip install h5py==2.10.0 --force-reinstall
conda deactivate
