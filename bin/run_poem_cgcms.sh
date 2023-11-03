#!/bin/bash
# Predict operons with short version of POEM
# run_poem_cgms.sh <path to input dir>

# get the scripts path
SCRIPT=`realpath -s $0`
SCRIPTPATH=`dirname $SCRIPT`

# set the path for the python
python=python

#######################################
# parse the args
######################################
while [ $# -gt 1 ]
do
    key="$1"
    case $key in
        -f|--fasta)
        temp="$2"
        shift # past argument
        ;;
        -a|--assembly)
        asm="$2"
        shift # past argument
        ;;
        -p|--predict)
        gpd="$2"
        shift # past argument
        ;;
        -l|--fasta)
        lr="$2"
        shift # past argument
        ;;
        *)
            # unknown option
        ;;
    esac
    shift # past argument or value
done

fasta=$temp/input.fsa
$python $SCRIPTPATH/../lib/prod2gmk.py $fasta\_prod_aa.fsa > $fasta\_gmk_aa.fsa
$python $SCRIPTPATH/../lib/reid.py $fasta\_gmk_aa.fsa > $fasta\_aa.fsa
$python $SCRIPTPATH/../lib/to_list.py $fasta\_aa.fsa > $fasta\.locus
$python $SCRIPTPATH/../lib/predict_operon.py predict $fasta $fasta\.locus $SCRIPTPATH/../config/Operon_Predictor/model.hdf5 > $fasta\.adjacency


