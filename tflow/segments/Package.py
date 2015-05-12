#TFLOW Segment: Copy and Zip Final Sequence File Result.
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 05/12/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import subprocess
import os.path
import sys

if __name__ == "__main__" and __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../'))
    import tflow.segments
    __package__ = "tflow.segments"

from .parser_class import OutputParser
from ..util import (print_exit, write_file, read_file)

JOB_TYPE = 'Package'
PROGRAM_URL = None
SEGMENT_FOR_VERSION = '0.9'
COMMAND_LIST = ['gzip']
COMMAND = ' '.join(COMMAND_LIST)
TEST_COMMAND = '-h'
OUT_FILE = JOB_TYPE + '.out'
MILESTONES = ['Packaging Final Sequence Output',
              'Sequence File Packaging Complete',
              ]
TERMINAL_FLAGS = []
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Not Found']

DEFAULT_SETTINGS = {'packaged_file_name':'Final_Assembly.fa',
                    #TFLOW Settings
                    'command_list':COMMAND_LIST,
                    'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    'write_command':True,
                    }
                    
REQUIRED_SETTINGS = ['out_file', 'command_list', 'write_command', 'working_directory',
                      'packaged_file_name']

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

def read(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.read_or_notify()

def test(options, silent=False):
    try:
        output = subprocess.check_output(options['command_list'] + [options['test_command']])
        if silent:
            return True
        else:
            print ' -- %s Found!' % JOB_TYPE
    
    except OSError as error:
        if silent:
            return False
        print '%s Cannot Be Found With Shell Command: "%s"' % (JOB_TYPE, options['command'])
        output = 'Error Number: %s\nError Text:\n%s' % (str(error.errno), error.strerror)
    return output

def analyze(options):
    print 'No Analysis Applicable to %s Segment.' % JOB_TYPE
    return ''

def run(options):
    if __name__ != '__main__' and options['is_pipe']:
        out_file = open(options['out_file'], 'w')
        terminal_out, terminal_error = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = out_file, out_file

    print 'Packaging Final Sequence Output'
    print ''

    #Ensure Required Settings in Options
    for required_option in REQUIRED_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for %s not given.' % (required_option, JOB_TYPE))

    #Ensure A Type of Input File is Given
    if not any(x in options for x in ['absolute_sequence_file',
                                      'rel_sequence_file',
                                      'result_name_file']):
        print_exit('Either absolute_sequence_file, rel_sequence_file, or'
                   + ' result_name_file paramater required.')

    #Ensure Working Directory Exists
    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    #Assign Correct Input Sequence File Name
    if 'absolute_sequence_file' in options:
        full_input_file = options['absolute_sequence_file']

    elif 'rel_sequence_file' in options:
        full_input_file = os.path.join(options['project_directory'],
                                       options['rel_sequence_file'])

    elif 'result_name_file' in options:
        full_result_name_file = os.path.join(options['project_directory'],
                                             options['result_name_file'])
        if os.path.isfile(full_result_name_file):
            print ('Reading Result Sequence File Name from Provided '
                   + 'File: %s' % full_result_name_file )
        else:
            print_exit('Provided File: %s Containing' % full_result_name_file
                       + ' Result Sequence File Name Not Found.')

        full_input_file = read_file(full_result_name_file)

        print 'Read Result Sequence File Name: %s' % full_input_file

    input_file = os.path.basename(full_input_file)

    #Check that Input Sequence File Exists
    if os.path.isfile(full_input_file):
        print 'Sequence Result File Found:'
        print '    ' + full_input_file
        print ''
    else:
        print_exit('Input Sequence File: %s Not Found.' % full_input_file)

    full_output_file = os.path.join(options['working_directory'], options['packaged_file_name'])

    if os.path.relpath(full_input_file, full_output_file) != '.':
        print 'Copying Final Result Sequence File:'
        print '    ' + full_input_file
        print 'To Location:'
        print '    ' + full_output_file
        print ''

        input_file_object = open(full_input_file, 'r')
        output_file_object = open(full_output_file, 'w')
        for line in input_file_object:
            output_file_object.write(line)
        input_file_object.close()
        output_file_object.close()

    print 'Zipping Final Result Sequence File:'
    print '    ' + full_output_file
    print ''

    zip_file_name = full_output_file + '.gz'
    command_list = list(options['command_list']) + ['-c', full_output_file]
    equivalent_command_list = list(command_list) + ['>', zip_file_name]
    command = ' '.join(equivalent_command_list)

    if options['write_command']:
        command_file = os.path.join(options['project_directory'],
                                    options['job_type'] + '.auto.sh')
        write_file(command_file, '#!/bin/sh\n' +command)

    print 'Running Command:\n    ' + command
    print ''

    sys.stdout.flush()
    zip_file = open(zip_file_name, 'w')
    process = subprocess.Popen(command_list, stdout=zip_file, stderr=sys.stderr, 
                                cwd=options['working_directory'])
    process.wait()
    zip_file.close()
    sys.stdout.flush()
                               
    if os.path.isfile(zip_file_name):
        print 'Zipped Sequence File: %s Created Successfully.' % zip_file_name
        print ''
    else:                              
        print_exit('Expected Output Zipped Sequence File: %s Cannot Be Found.' % zip_file_name)
                                   
    print 'Sequence File Packaging Complete'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout, sys.stderr = terminal_out, terminal_error
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
