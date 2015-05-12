#!/bin/bash
#Convenience Wrapper for Main TFLOW Component: tflow/manifold.py
#Usage: "tflow.sh RUN_MODE [--options]"
#For Advanced Usage: "tflow.sh -h"
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

TFLOW_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
#export PYTHONPATH=$TFLOW_DIR:$PYTHONPATH
python2.7 $TFLOW_DIR/tflow/manifold.py "$@"

