#TFLOW Segment: Parse provided reads files into lists, when possible.
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import subprocess
import os.path
import sys

if __name__ == "__main__" and __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../'))
    import tflow.segments
    __package__ = "tflow.segments"

from .parser_class import OutputParser
from ..util import (print_exit, write_file, write_file_list, ensure_FASTQ_GZ, 
                    ensure_FASTA_GZ, ensure_list)
from .. import util

JOB_TYPE = 'Make_Read_Lists'
PROGRAM_URL = None
SEGMENT_FOR_VERSION = '0.9'
#COMMAND_LIST = ['python2.7', os.path.realpath(__file__)]
#COMMAND = ' '.join(COMMAND_LIST)
#TEST_COMMAND = '-h'
OUT_FILE = JOB_TYPE + '.out'
MILESTONES = ['Preparation Complete']
TERMINAL_FLAGS = []
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Exception: ERROR',
                 'Not Found']

DEFAULT_SETTINGS = {'is_paired_reads':True,
                    'read_type':'fq',
                    'raw_single_reads_list':'raw_unpaired_reads_list',
                    'raw_left_reads_list':'raw_left_reads_list',
                    'raw_right_reads_list':'raw_right_reads_list',
                    'left_read_indicator':'_R1',
                    'right_read_indicator':'_R2',
                    #TFLOW Settings
                    #'command':COMMAND,
                    #'command_list':COMMAND_LIST,
                    #'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    'write_command':True,
                    }
                    
#REQUIRED_SETTINGS = ['command_list', 'is_paired_reads', 'read_type', 'write_command']
REQUIRED_SETTINGS = ['is_paired_reads', 'read_type']

class Parser(OutputParser):
    def set_local_defaults(self):
        self.milestones = MILESTONES
        self.terminal_flags = TERMINAL_FLAGS
        self.failure_flags = FAILURE_FLAGS
        self.job_type = JOB_TYPE

def check_done(options):
    parser = Parser()
    parser.out_file = options['out_file']
    failure_exit = (options['mode'] in ['run', 'track'])
    return parser.check_completion(failure_exit)

def track(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.track()

def analyze(options):
    print '    Analysis Not Applicable '

def read(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.read_or_notify()

def stop(options):
    print '    Job Stopping Not Applicable'

def clean(options):
    remove_outfile = (options['mode'] == 'reset')
    files = [options['raw_single_reads_list'], options['raw_left_reads_list'],
             options['raw_right_reads_list']]
    util.clean_TFLOW_auto_files(options['job_type'], options['project_directory'],
                                options['working_directory'], remove_outfile=remove_outfile,
                                confirm=options['confirm'], files=files)

def test(options, silent=False):
    if silent:
        return True
    else:
        print ' -- %s Found!' % JOB_TYPE
        output = 'File Location: %s' % os.path.realpath(__file__)
        return output

    #Individual Executable Operation Currently Disabled.
    #try:
    #    output = subprocess.check_output(options['command_list'] + [options['test_command']])
    #    if silent:
    #        return True
    #    else:
    #        print ' -- %s Found!' % JOB_TYPE

    #except OSError as error:
    #    if silent:
    #        return False
    #    print '%s Cannot Be Found With Shell Command: "%s"' % (JOB_TYPE, options['command'])
    #    output = 'Error Number: %s\nError Text:\n%s' % (str(error.errno), error.strerror)
    #return output

def run(options):
    if __name__ != '__main__' and options['is_pipe']:
        out_file = open(options['out_file'], 'w')
        terminal_out = sys.stdout
        terminal_error = sys.stderr
        sys.stdout = out_file
        sys.stderr = out_file

    for required_option in REQUIRED_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for Make_Read_Lists not given.' % required_option)
        
    if 'raw_reads' in options:
        reads = ensure_list(options['raw_reads'])
        print 'Parsing Input Reads:'
        for read in reads:
            print '  -- ' + read

        if options['is_paired_reads']:
            left_reads = []
            right_reads = []
            for read in reads:
                if options['left_read_indicator'] in read:
                    left_reads.append(read)
                elif options['right_read_indicator'] in read:
                    right_reads.append(read)
                else:
                    print_exit('Pairing of Read %s Cannot Be Established' % read
                               + 'Using Pairing Indicators '
                               + 'Left: "%s" ' % options['left_read_indicator']
                               + ' and Right: %s.\n' % options['right_read_indicator']
                               + 'Please set indicators.')

    elif all(x in options for x in ['raw_left_reads', 'raw_right_reads']):
        left_reads = ensure_list(options['raw_left_reads'])
        right_reads = ensure_list(options['raw_right_reads'])
        reads = left_reads + right_reads
        print 'Parsing Paired Input Reads:'
        print 'Left Reads:'
        for read in left_reads:
            print '  -- ' + read
        print ''
        print 'Right Reads:'
        for read in right_reads:
            print '  -- ' + read
        print ''
        
    else:
        if options['is_paired_reads']:
            print_exit('Required Options: left_reads, right_reads, for input'
                       + ' of paired end reads not found.')
        else:
            print_exit('Required Option: reads, for input'
                       + ' of paired end reads not found.')

    #Individual Executable Operation Currently Disabled.
    #if options['is_paired_reads']:
    #    command = (COMMAND + ' --raw_left_reads %s' % ' '.join(left_reads)
    #               + ' --raw_right_reads %s ' % ' '.join(right_reads)
    #               + '--is_paired_reads True') 
    #else:
    #    command = (COMMAND + ' --reads %s' % ' '.join(reads)
    #               + '--is_paired_reads False') 
    #
    #if options['write_command']:
    #    command_file = os.path.join(options['project_directory'], 
    #                                options['job_type'] + '.auto.sh')
    #    write_file(command_file, '#!/bin/sh\n' + command)
    #print 'Running Command:\n    ' + command
    #print ''

    #Ensure reads exist and are of correct type.
    for read in reads:
        full_read = os.path.join(options['working_directory'], read)

        if options['read_type'] == 'fq':
            ensure_FASTQ_GZ(full_read)
        elif options['read_type'] == 'fa':
            ensure_FASTA_GZ(full_read)
        else:
            print_exit('Provided read_type value %s not "fq" or "fa"' % options['read_type'])

        print '  Read: %s Found.' % read

    print ''

    if options['is_paired_reads']:
        for required_option in ['raw_left_reads_list', 'raw_right_reads_list']:
            if required_option not in options:
                print_exit('Required Option: %s for Make_Read_Lists not given.' % required_option)

        if len(left_reads) != len(right_reads):
            print_exit('Number of Left Reads: %i' % len(left_reads)
                       + 'Does Not Equal Number of Right Reads: %i' % len(right_reads))

        print 'Identified Paired End Reads:'
        write_file_list(options['raw_left_reads_list'], left_reads)
        print '  Left read file names written to file: %s' % options['raw_left_reads_list']
        for read in left_reads:
            print '  -- ' + read
        print ''
        write_file_list(options['raw_right_reads_list'], right_reads)
        print '  Right read file names written to file: %s' % options['raw_right_reads_list']
        for read in right_reads:
            print '  -- ' + read

    #Unpaired Reads
    else:
        if 'raw_single_reads_list' not in options:
            print_exit('Required Option: raw_single)reads_list for Make_Read_Lists not given.')

        print 'Using Unpaired Reads:'
        write_file_list(options['raw_single_reads_list'], reads)
        print '  Unpaired read names written to file: %s' % options['raw_single_reads_list']
        for read in reads:
            print '  -- ' + read

    print ''
    print 'Read File List Preparation Complete.'                        

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout = terminal_out
        sys.stderr = terminal_error
        out_file.close()

#Individual Executable Operation Currently Disabled.
#if __name__ == '__main__':
#    from ..manifold import get_settings, print_settings
#    options = get_settings()
#    if 'project_directory' not in options:
#        options['project_directory'] = os.getcwd()
#    if 'working_directory' not in options:
#        options['working_directory'] = os.getcwd()
#    #print_settings(options)
#
#    print ''
#    print 'Running %s Job...' % JOB_TYPE
#    print ''
#
#    run(options)
#
#    print ''
#    print 'Complete'
#    print ''
