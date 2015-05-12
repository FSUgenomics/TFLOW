#!/bin/bash
#Convenience Function to Test TFLOW Setup
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

TFLOW_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#export PYTHONPATH=$TFLOW_DIR:$PYTHONPATH
$TFLOW_DIR/tflow.sh test -t Test_Pipe 
#$TFLOW_DIR/tflow.sh test -t Test_Pipe --print_test_output=True
