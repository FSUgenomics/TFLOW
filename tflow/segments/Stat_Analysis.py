#TFLOW Segment: Perform Statistical Analysis on FASTA File Sequences
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 05/12/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import os.path
import sys
import subprocess
import shutil

if __name__ == "__main__" or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../'))
    import tflow.segments
    __package__ = "tflow.segments"

from .parser_class import OutputParser
from ..fasta import check_N50_in_place
from ..util import print_exit, read_file, write_file, write_report

JOB_TYPE = 'Stat_Analysis'
PROGRAM_URL = None
SEGMENT_FOR_VERSION = '0.9'
#COMMAND_LIST = ['python2.7', os.path.realpath(__file__)]
#COMMAND = ' '.join(COMMAND_LIST)
#TEST_COMMAND = '-h'
OUT_FILE = JOB_TYPE + '.out'
MILESTONES = ['Performing Statistical Analysis',
              'Statistical Analysis Complete']
TERMINAL_FLAGS = []
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Not Found']
DEFAULT_SETTINGS = {'copy_input_file':False,
                    #TFLOW Settings
                    #'command':COMMAND,
                    #'command_list':COMMAND_LIST,
                    #'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    'write_report':True,
                    #'write_command':True,
                    }

REQUIRED_SETTINGS = ['working_directory', 'copy_input_file', 'write_report']
REQUIRED_ANALYSIS_SETTINGS = REQUIRED_SETTINGS

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
    #
    #except OSError as error:
    #    if silent:
    #        return False
    #    print '%s Cannot Be Found With Shell Command: "%s"' % (JOB_TYPE, options['command'])
    #    output = 'Error Number: %s\nError Text:\n%s' % (str(error.errno), error.strerror)
    #return output

def run(options):
    if __name__ != '__main__' and options['is_pipe']:
        out_file = open(options['out_file'], 'w')
        terminal_out, terminal_error = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = out_file, out_file

    print 'Performing Statistical Analysis'
    print ''
    analyze(options)

    print 'Statistical Analysis Complete'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout, sys.stderr = terminal_out, terminal_error
        out_file.close()

#Analyze Results of Sequence Comparison
def analyze(options):
    analysis = ''

    #Ensure Required Settings in Options
    for required_option in REQUIRED_ANALYSIS_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for %s not given.' % (required_option, JOB_TYPE))

    #Ensure Working Directory Exists
    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    #Ensure A Type of Input File is Given
    if not any(x in options for x in ['absolute_input_analysis_file',
                                      'rel_input_analysis_file',
                                      'result_name_file']):
        print_exit('Either absolute_input_analysis_file, rel_input_analysis_file, or'
                   + ' result_name_file paramater required.')

    #Assign Correct Input File Name
    if 'absolute_input_analysis_file' in options:
        full_input_file = options['absolute_input_analysis_file']

    elif 'rel_input_analysis_file' in options:
        full_input_file = os.path.join(options['project_directory'],
                                       options['rel_input_analysis_file'])

    elif 'result_name_file' in options:
        full_result_name_file = os.path.join(options['project_directory'],
                                             options['result_name_file'])
        if os.path.isfile(full_result_name_file):
            analysis += ('Reading Result Sequence File Name from Provided '
                         + 'File: %s\n' % full_result_name_file )
        else:
            print_exit('Provided File: %s Containing' % full_result_name_file
                       + ' Result Sequence File Name Not Found.')

        full_input_file = read_file(full_result_file_name)

    if not os.path.isfile(full_input_file):
        print_exit('Cannot Find Read Result Sequence File: %s' % full_input_file)

    input_file = os.path.basename(full_input_file)
    analysis += 'Performing Statistical Analysis on File:\n'
    analysis += '    %s\n' % full_input_file

    results = check_N50_in_place(full_input_file, fail_exit=False, return_report_dict=True)
    if isinstance(results, tuple):
        report_analysis, report_dict = results
        analysis += report_analysis
    else:
        analysis += results
        report_dict = None

    if options['write_report'] and report_dict:
        report_file = os.path.join(options['working_directory'],
                                   JOB_TYPE + '.report')
        write_report(report_file, report_dict)

    return analysis

            
    
    
    

