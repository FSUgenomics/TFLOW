#TFLOW Segment: Assemble Sequences with CAP3 Assembler
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow


#!/usr/bin/env python
#Cap3 Sequencing Tool Parameters
#Dan Stribling
#Version 1.0, 09/21/2014

import os.path
import sys
import subprocess
import shutil

if __name__ == "__main__" or __package__ is None:
    import tflow.segments
    __package__ = "tflow.segments"

from .. import local_settings
from .parser_class import OutputParser
from ..util import print_exit, print_warning, write_file
from ..fasta import label_sequences, check_N50_in_place

if hasattr(local_settings, 'CAP3_LOCATION'):
    CAP3_LOCATION = local_settings.CAP3_LOCATION
else:
    CAP3_LOCATION = ''

if hasattr(local_settings, 'CAP3_EXEC'):
    CAP3_EXEC = local_settings.CAP3_EXEC
else:
    CAP3_EXEC = os.path.join(CAP3_LOCATION,'cap3')


JOB_TYPE = 'CAP3'
PROGRAM_URL = 'http://doua.prabi.fr/software/cap3'
SEGMENT_FOR_VERSION = 'N/A'
COMMAND = CAP3_EXEC
COMMAND_LIST = [COMMAND]
TEST_COMMAND = None
OUT_FILE = 'CAP3.out'
MILESTONES = ['CAP3 Job Done']
TERMINAL_FLAGS = []
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Not Found']
DEFAULT_SETTINGS = {'working_directory':'CAP3',
                    'connections_file':'connections.out',
                    #TFLOW Trinity Settings
                    'command':COMMAND,
                    'command_list':COMMAND_LIST,
                    'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    'write_report':True,
                    'write_command':True,
                    }

REQUIRED_SETTINGS = ['command_list', 'working_directory', 'write_report', 'write_command']

REQUIRED_ANALYSIS_SETTINGS = ['working_directory', 'write_report']


class Parser(OutputParser):
    def set_local_defaults(self):
        self.milestones = MILESTONES
        self.terminal_flags = TERMINAL_FLAGS
        self.failure_flags = FAILURE_FLAGS
        self.job_type = JOB_TYPE

def check_done(options):
    parser = Parser()
    parser.out_file = options['out_file']
    return parser.check_completion()

def track(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.track()

def analyze(options):
    for required_option in REQUIRED_ANALYSIS_SETTINGS:
        if required_option not in options:
            print_warning('Required Option: %s for ' % required_option
                          + '%s Analysis not given.' % JOB_TYPE)
            return ''

    if not any(x in options for x in ['absolute_input_file', 'relative_input_file']):
        print_warning('Either absolute_input_file or relative_input_file paramater required.')
        return ''

    if not os.path.isdir(options['working_directory']):
        print_warning('Working Directory: %s Not Found.' % options['working_directory'])
        return ''

    if 'absolute_input_file' in options:
        full_input_file = options['absolute_input_file']
        input_file = os.path.basename(full_input_file)

    elif 'relative_input_file' in options:
        full_input_file = os.path.join(options['project_directory'], 
                                       options['relative_input_file'])
        input_file = os.path.basename(options['relative_input_file'])


    expected_full_output_file = os.path.join(options['working_directory'],
                                             input_file + '.cap.combined')
    if not os.path.isfile(expected_full_output_file):
        print_warning('Expected Output File: %s Not Found.' % expected_full_output_file)
        return ''

    results = check_N50_in_place(expected_full_output_file, fail_exit=False, return_report=True)
    if isinstance(results, tuple):
        analysis, report = results
    else:
        analysis = results
        report = None

    if options['write_report'] and report:
        report_file = os.path.join(options['working_directory'],
                                   JOB_TYPE + '.report')
        write_file(report_file, report)

    return analysis



def read(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.read_or_notify()

def test(options, silent=False):
    try:
        process = subprocess.Popen(options['command_list'], stdout=subprocess.PIPE, 
                         stderr=subprocess.STDOUT)
        process.wait()
        output, error = process.communicate()
        print ' -- %s Found!' % JOB_TYPE

    except OSError as error:
        if silent:
            return False
        print '%s Cannot Be Found With Shell Command: "%s"' % (JOB_TYPE, options['command'])
        if PROGRAM_URL:
            print 'If Not Installed, %s Can be Downloaded From:\n%s' % (JOB_TYPE, PROGRAM_URL)
        output = 'Error Number: %s\nError Text:\n%s' % (str(error.errno), error.strerror)

    return output

def run(options):
    if __name__ != '__main__' and options['is_pipe']:
        out_file_stream = open(options['out_file'], 'w')
        terminal_out, terminal_error = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = out_file_stream, out_file_stream

    for required_option in REQUIRED_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for %s not given.' % (required_option, JOB_TYPE))

    if not any(x in options for x in ['absolute_input_file', 'relative_input_file']):
        print_exit('Either absolute_input_file or relative_input_file paramater required.')

    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    if 'absolute_input_file' in options:
        full_input_file = options['absolute_input_file']
        input_file = os.path.basename(full_input_file)

    elif 'relative_input_file' in options:
        full_input_file = os.path.join(options['project_directory'], 
                                       options['relative_input_file'])
        input_file = os.path.basename(options['relative_input_file'])

    if not os.path.isfile(full_input_file):
        print_exit('Input File: %s Not Found.' % full_input_file)

    print 'Copying Input File: %s to Working Directory: %s' % (input_file, 
                                                               options['working_directory']) 
    working_input_file = os.path.join(options['working_directory'], input_file)
    shutil.copyfile(full_input_file, working_input_file)

    if not os.path.isfile(working_input_file):
        print_exit('Copying of File: %s to Name: %s Unsuccesful.' % (full_input_file, 
                                                                     working_input_file))

    command_list = options['command_list'][:]
    command_list.append(working_input_file)
    for possible_option in [('-' + x) for x in ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
                                                'k', 'm', 'n', 'o', 'p', 'r', 's', 't', 'u', 'v',
                                                'w', 'x' 'y', 'z']]:
        if possible_option in options:
            command_list += [possible_option, options[possible_option]]

    command = ' '.join(command_list)

    if options['write_command']:
        command_file = os.path.join(options['project_directory'],
                                    options['job_type'] + '.auto.sh')
        write_file(command_file, '#!/bin/sh\n' +command)

    print ''
    print 'Running Command:\n    ' + command

    sys.stdout.flush()

    if 'connections_file' in options:
        full_connections_file = os.path.join(options['working_directory'], 
                                             options['connections_file'])
        process_out = open(full_connections_file, 'w')
        print 'Writing CAP3 Connection Information to File: %s' % full_connections_file
    else:
        process_out = sys.stdout

    try:
        process = subprocess.Popen(command_list, stdout=process_out, stderr=sys.stderr,
                                   cwd=options['project_directory'])
        process.wait()
        sys.stdout.flush()

        if process_out != sys.stdout:
            process_out.close()

    except KeyboardInterrupt:
        if __name__ != '__main__' and options['is_pipe']:
            sys.stdout, sys.stderr = terminal_out, terminal_error
            out_file_stream.close()

        print 'Killing %s Process.' % JOB_TYPE
        process.kill()
        raise

    expected_output_contigs = os.path.join(options['working_directory'], 
                                           input_file + '.cap.contigs')
    expected_output_singlets = os.path.join(options['working_directory'], 
                                              input_file + '.cap.singlets')
    
    for expected_output in [expected_output_contigs, expected_output_singlets]:
        if not os.path.isfile(expected_output):
            print_exit('Expected Output %s Not Found!' % expected_output)

    print 'All CAP3 Outputs Found!'
    if 'label' in options and options['label']:
        print 'Labeling Singlets and Contigs from Output Files with Label: %s' % options['label']

        label_sequences(expected_output_contigs, options['label'])
        label_sequences(expected_output_singlets, options['label'] + '_singlets')
        expected_output_contigs += '.labeled'
        expected_output_singlets += '.labeled'

    concatenated_output = os.path.join(options['working_directory'],
                                       input_file + '.cap.combined')
    print 'Creating Concatenated Output CAP3 Sequence File: %s' % concatenated_output
    cat_file = open(concatenated_output, 'w')
    for expected_output in [expected_output_contigs, expected_output_singlets]:
        expected_output_file = open(expected_output, 'r')
        for line in expected_output_file:
            cat_file.write(line)
        expected_output_file.close()
    cat_file.close()

    analyze(options)

    print ''
    print 'CAP3 Job Done'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout = terminal_out
        sys.stderr = terminal_error
        out_file_stream.close()