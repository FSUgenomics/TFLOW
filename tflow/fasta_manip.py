#!/usr/bin/env python2.7
#TFLOW Utility: Manipulate or Analyze FASTA File.
#Usage: "fasta_manip.py MODE sequence_file.fa"
#For Full Usage: "fasta_manip.py -h"
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import argparse
import sys
import os.path

if __name__ == "__main__" and __package__ is None:
    import tflow
    __package__ = "tflow"

from .fasta import check_FASTA, check_N50, check_N50_in_place

PRINT_MODES = {'check':'check', 'N50':'N50', 'IP_N50':'IP_N50'}
FLEXIBLE_MODES = {'verify':'check', 'n50':'N50', 'ip_n50':'IP_N50', 'IP_n50':'IP_N50', 
                  'ip_N50':'IP_N50'}
ALL_MODES = PRINT_MODES.copy()
ALL_MODES.update(FLEXIBLE_MODES)

def flexible_mode(mode):
    if mode in ALL_MODES:
        return ALL_MODES[mode]
    return mode

def parse_count_sequences_args():
    parser = argparse.ArgumentParser(prog='fasta_manip.py',
                                     description='Analyze and Manipulate FASTA Sequence File')

    parser.add_argument('mode', action='store', help='Run Mode', type=flexible_mode, 
                        choices=PRINT_MODES.keys())
    parser.add_argument('sequence_file', action='store', 
                        help='Input Sequence File', metavar='SEQUENCE_FILE')

    return vars(parser.parse_args())

if __name__ == '__main__':

    options = parse_count_sequences_args()
    sequence_file = options['sequence_file']
    mode = options['mode']

    if not os.path.isfile(sequence_file):
        print 'File: %s Not Found.\nExiting...' % sequence_file
        sys.exit(1) 

    if mode == 'check':
        print ''
        print ' --- Checking FASTA file %s for Formatting Issues ---' % sequence_file
        print ''
        check_FASTA(sequence_file)

    elif mode == 'N50':
        print ''
        print '  --- Going to Find N50 Length of Fasta File: %s ---' % sequence_file
        print ''
        check_N50(sequence_file)

    elif mode == 'IP_N50':
        print ''
        print '  --- Going to Find N50 Length of Fasta File: %s In Place ---' % sequence_file
        print ''
        check_N50_in_place(sequence_file)

    print ''
    print ' --- All Done! --- '
    print ''
        
