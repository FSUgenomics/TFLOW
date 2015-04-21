#!/bin/bash
#Convenience Wrapper for TFLOW Utility: tflow/count_sequences.py
#Count sequences in one or more .fasta[.gz], .fastq[.gz], .fna etc... files.
#Usage: "count_sequences.sh [sequence_file_1.fa] [sequence_file....]"
#For Full Usage: "count_sequences.sh -h"
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

TFLOW_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
export PYTHONPATH=$TFLOW_DIR:$PYTHONPATH
python2.7 $TFLOW_DIR/tflow/count_sequences.py "$@"

