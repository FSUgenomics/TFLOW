#!/bin/bash
#Convenience Wrapper for TFLOW Utility: tflow/label_sequences.py
#Add prefix label to sequence headers in a FASTA file.
#Usage: "label_sequences.sh input_file.fa prefix_label [output_file.fa]"
#For Full Usage: "label_sequences.sh -h"
#
#Dan Stribling
#FSU Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

TFLOW_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#export PYTHONPATH=$TFLOW_DIR:$PYTHONPATH
python2.7 $TFLOW_DIR/tflow/label_sequences.py "$@"

