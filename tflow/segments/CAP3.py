#TFLOW Segment: Assemble Sequences with CAP3 Assembler
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import os.path
import sys
import subprocess
import shutil
from time import sleep

if __name__ == "__main__" or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../'))
    import tflow.segments
    __package__ = "tflow.segments"

from .. import local_settings
from .parser_class import OutputParser
from ..util import (print_exit, print_warning, write_file, write_report, read_file, 
                    delete_pid_file, stop_TFLOW_process)
from .. import util
from ..fasta import label_sequences, check_N50_in_place, check_FASTA

if hasattr(local_settings, 'CAP3_LOCATION'):
    CAP3_LOCATION = local_settings.CAP3_LOCATION
else:
    CAP3_LOCATION = ''

if hasattr(local_settings, 'CAP3_EXEC'):
    CAP3_EXEC = local_settings.CAP3_EXEC
else:
    CAP3_EXEC = os.path.join(CAP3_LOCATION,'cap3')


JOB_TYPE = 'CAP3'
PROGRAM_URL = 'http://seq.cs.iastate.edu/cap3.html'
SEGMENT_FOR_VERSION = 'N/A'
COMMAND = CAP3_EXEC
COMMAND_LIST = [COMMAND]
TEST_COMMAND = None
OUT_FILE = 'CAP3.out'
MILESTONES = ['CAP3 Job Done']
TERMINAL_FLAGS = []
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Exception: ERROR',
                 'Not Found']
DEFAULT_SETTINGS = {'working_directory':'CAP3',
                    'connections_file':'connections.out',
                    'combined_input_name':'combined_input.fa',
                    #TFLOW Trinity Settings
                    'command':COMMAND,
                    'command_list':COMMAND_LIST,
                    'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    'write_report':True,
                    'write_command':True,
                    'write_result_name':False,
                    'write_pid':True,
                    }

REQUIRED_SETTINGS = ['command_list', 'working_directory', 'write_report', 'write_command', 
                     'write_pid', 'combined_input_name', 'write_result_name']

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
    failure_exit = (options['mode'] in ['run', 'track'])
    return parser.check_completion(failure_exit)

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

    #Ensure that Input File Options Given
    if any(x in options for x in ['absolute_input_file', 'relative_input_file',
                                  'absolute_input_files', 'relative_input_files']):
        pass

    elif options['write_result_name']:
        report_name_file_name = os.path.join(options['working_directory'],
                                             JOB_TYPE + '.auto.result_name')
        if not os.path.isfile(report_name_file_name):
            print_warning('Either absolute_input_file, relative_input_file, absolute_input_files,'
                          + ' or relative_input_files paramater required.')
            return ''

        print ('Reading Result Sequence File Name from Provided '
               + 'File: %s' % report_name_file_name )

        full_input_file = read_file(report_name_file_name)
        if os.path.isfile(full_input_file):
            print 'Read Result Sequence File Name: %s' % full_input_file
            print 'File Found!'
            print ''
        else:
            print 'Cannot Find Read Result Sequence File: %s' % full_input_file
            return ''

    else:
        print_warning('Either absolute_input_file, relative_input_file, absolute_input_files,'
                      + ' or relative_input_files paramater required.')
        return ''
 
    if not os.path.isdir(options['working_directory']):
        print_warning('Working Directory: %s Not Found.' % options['working_directory'])
        return ''

    if 'absolute_input_file' in options:
        full_input_file = options['absolute_input_file']

    elif 'relative_input_file' in options:
        full_input_file = os.path.join(options['project_directory'], 
                                       options['relative_input_file'])

    elif any(x in options for x in ['absolute_input_files', 'relative_input_files']):
        full_input_file = os.path.join(options['working_directory'],
                                       options['combined_input_name'])

    input_file = os.path.basename(full_input_file)
    expected_full_output_file = os.path.join(options['working_directory'],
                                             input_file + '.cap.combined')
    if not os.path.isfile(expected_full_output_file):
        print_warning('Expected Output File: %s Not Found.' % expected_full_output_file)
        return ''

    results = check_N50_in_place(expected_full_output_file, fail_exit=False, 
                                 return_report_dict=True)
    if isinstance(results, tuple):
        analysis, report_dict = results
    else:
        analysis = results
        report_dict = None

    if options['write_report'] and report_dict:
        report_file = os.path.join(options['working_directory'],
                                   JOB_TYPE + '.report')
        write_report(report_file, report_dict)

    return analysis



def read(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.read_or_notify()

def stop(options):
    job_pid_file = os.path.join(options['working_directory'],
                                JOB_TYPE + '.auto.pid')
    stop_TFLOW_process(job_pid_file, JOB_TYPE)

def clean(options):
    out_files = ['connections.out', 'combined_input.fa']
    remove_outfile = (options['mode'] == 'reset')
    util.clean_TFLOW_auto_files(options['job_type'], options['project_directory'], 
                                options['working_directory'], remove_outfile=remove_outfile, 
                                confirm=options['confirm'], out_files=out_files)

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

    if not any(x in options for x in ['absolute_input_file', 'relative_input_file',
                                      'absolute_input_files', 'relative_input_files']):
        print_exit('Either absolute_input_file, relative_input_file, absolute_input_files,'
                   + ' or relative_input_files paramater required.')

    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    if 'absolute_input_file' in options:
        full_input_file = options['absolute_input_file']

    elif 'relative_input_file' in options:
        full_input_file = os.path.join(options['project_directory'], 
                                       options['relative_input_file'])

    elif any(x in options for x in ['absolute_input_files', 'relative_input_files']):
        if 'absolute_input_files' in options:
            input_files = list(options['absolute_input_files'])
        elif 'relative_input_files' in options:
            print 'Relative Input Files Provided:'
            for input_file in options['relative_input_files']:
                print '   ', input_file
            print ''
            input_files = [os.path.join(options['project_directory'], x) 
                           for x in options['relative_input_files']]

        print 'Creating Combined Input File %s' % options['combined_input_name']
        print '    From Input Files:'
        for input_file in input_files:
            print '   ', input_file      
        print ''

        combined_seq_file_name = os.path.join(options['working_directory'],
                                              options['combined_input_name'])
        with open(combined_seq_file_name, 'w') as combined_seq_file:
            for input_file_name in input_files:
                print 'Adding Sequences from File: %s' % input_file_name
                sys.stdout.flush()
                input_file = open(input_file_name, 'r')
                for line in input_file:
                    combined_seq_file.write(line)
                combined_seq_file.flush()
                input_file.close()

        print 'Creation of Combined Input File: %s Completed.' % combined_seq_file_name
        print ''
        #print 'Skipping Verification of Integrity of Combined FASTA File:'
        print 'Verifying Integrity of Combined FASTA File:'
        #from time import sleep
        #sleep(3)
        check_FASTA(combined_seq_file_name)
        print 'Verification of FASTA Complete.'
        print ''
        
        full_input_file = combined_seq_file_name

    input_file = os.path.basename(full_input_file)

    if not os.path.isfile(full_input_file):
        print_exit('Input File: %s Not Found.' % full_input_file)

    working_input_file = os.path.join(options['working_directory'], input_file)
    if not os.path.isfile(working_input_file):
        print 'Copying Input File: %s to Working Directory: %s' % (input_file, 
                                                                   options['working_directory']) 
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
        if options['write_pid']:
            pid_file_name = os.path.join(options['working_directory'],
                                         options['job_type'] + '.auto.pid')
            write_file(pid_file_name, str(process.pid))
        process.wait()
        sys.stdout.flush()
        print 'Ensuring Wrapup...'
        sleep(60)

        if options['write_pid']:
            delete_pid_file(pid_file_name)

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
        label_sequences(expected_output_singlets, options['label'])
        expected_output_contigs += '.labeled'
        expected_output_singlets += '.labeled'

    concatenated_output = os.path.join(options['working_directory'],
                                       input_file + '.cap.combined')
    print 'Creating Concatenated Output CAP3 Sequence File: %s' % concatenated_output
    with open(concatenated_output, 'w') as cat_file:
        for expected_output in [expected_output_contigs, expected_output_singlets]:
            expected_output_file = open(expected_output, 'r')
            for line in expected_output_file:
                cat_file.write(line)
            expected_output_file.close()

    if options['write_result_name']:
        report_name_file_name = os.path.join(options['working_directory'],
                                             JOB_TYPE + '.auto.result_name')
        write_file(report_name_file_name, concatenated_output)


    analyze(options)

    print ''
    print 'CAP3 Job Done'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout = terminal_out
        sys.stderr = terminal_error
        out_file_stream.close()
