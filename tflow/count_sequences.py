#!/usr/bin/env python2.7
#TFLOW Utility: Count sequences in one or more .fasta[.gz], .fastq[.gz], .fna etc... files.
#Usage: "count_sequences.py [sequence_file_1.fa] [sequence_file....]"
#For Full Usage: "count_sequences.py -h"
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import argparse
import os

if __name__ == "__main__" and __package__ is None:
    import sys
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
    import tflow
    __package__ = "tflow"

#from .util import count_sequences, is_sequence_file, print_exit
from . import util

def parse_count_sequences_args():
    parser = argparse.ArgumentParser(prog='count_sequences.py',
                                     description='Count Sequences in Sequence Files')
    parser.add_argument('contents', action='store', nargs='*', default=[], 
                        help='Input Sequence File(s)', metavar='SEQUENCE_FILE')
    return vars(parser.parse_args())

if __name__ == '__main__':
    options = parse_count_sequences_args()
    contents = options['contents']
    for item in contents:
        if '*' in item:
            util.print_exit(['Argument: ( %s ) contains invalid symbol "*".',
                             'Was file name input correctly?',
                             ''], 1)
                       
    if not contents:
        contents = []
        for file_name in os.listdir(os.getcwd()):
            if util.is_sequence_file(file_name):
                contents.append(file_name)

    if not contents:
        util.print_exit(['No Sequence Files Found.', ''], 1)

    util.count_sequences(contents)
