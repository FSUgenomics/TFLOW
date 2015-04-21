#!/bin/bash
#Convenience Wrapper for TFLOW Utility: tflow/fasta_manip.py
#Manipulate or Analyze FASTA File
#Usage: "fasta_manip.py MODE sequence_file.fa"
#For Full Usage: "fasta_manip.py -h"
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

TFLOW_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
export PYTHONPATH=$TFLOW_DIR:$PYTHONPATH
python2.7 $TFLOW_DIR/tflow/fasta_manip.py "$@"

