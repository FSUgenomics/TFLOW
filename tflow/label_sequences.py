#!/usr/bin/env python2.7
#TFLOW Utility: Add prefix label to sequence headers in a FASTA file.
#Usage: "label_sequences.py input_file.fa prefix_label [output_file.fa]"
#For Full Usage: "label_sequences.py -h"
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import argparse

if __name__ == "__main__" and __package__ is None:
    import tflow
    __package__ = "tflow"

from .util import ensure_FASTA
from .fasta import label_sequences

def parse_label_sequences_args():
    parser = argparse.ArgumentParser(prog='label_sequences.py',
                                     description='Add a Prefix Labels to FASTA sequences.')
    parser.add_argument('input_file_name', action='store', help='Input FASTA Sequence File', 
                        metavar='IN_FILE')
    parser.add_argument('label', action='store', help='Sequence Prefix Label',  metavar='LABEL')
    parser.add_argument('output_file_name', action='store', default=None,
                        nargs='?', help='Ouput File Name', metavar='OUT_FILE_NAME')
    return vars(parser.parse_args())


if __name__ == '__main__':
    
    options = parse_label_sequences_args()
    ensure_FASTA(options['input_file_name'])

    print ''
    label_sequences(options['input_file_name'], options['label'], options['output_file_name'])
    print ''
    print 'Labeling Complete'
    print ''
    
