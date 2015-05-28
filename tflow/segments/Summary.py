#TFLOW Segment: Summarize Reports from any Run TFLOW Segments.
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
from ..util import (print_exit, write_file, read_report, combine_report, write_file_list, 
                    ensure_FASTQ_GZ, ensure_FASTA_GZ, ensure_list)
from .. import util

JOB_TYPE = 'Summary'
PROGRAM_URL = None
SEGMENT_FOR_VERSION = '0.9'
#COMMAND_LIST = ['python2.7', os.path.realpath(__file__)]
#COMMAND = ' '.join(COMMAND_LIST)
#TEST_COMMAND = '-h'
OUT_FILE = JOB_TYPE + '.out'
MILESTONES = ['Creating Summary Report',
              'Summary Complete']
TERMINAL_FLAGS = []
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Exception: ERROR',
                 'Not Found']

DEFAULT_SETTINGS = {'write_report':True,
                    'write_csv_report':True,
                    #TFLOW Settings
                    #'command':COMMAND,
                    #'command_list':COMMAND_LIST,
                    #'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    #'write_command':True,
                    }
                    
#REQUIRED_SETTINGS = ['out_file', 'command_list', 'write_report', 'write_command', 
#                     'working_directory']
REQUIRED_SETTINGS = ['out_file', 'write_report', 'write_csv_report', 'working_directory']

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

def stop(options):
    print '    %s Job Stopping Not Applicable' % JOB_TYPE

def clean(options):
    remove_outfile = (options['mode'] == 'reset')
    util.clean_TFLOW_auto_files(options['job_type'], options['project_directory'],
                                options['working_directory'], remove_outfile=remove_outfile, 
                                confirm=options['confirm'])

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

def analyze(options):
    for required_option in REQUIRED_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for segment: ' % required_option
                       + '%s not given.' % JOB_TYPE )

    analysis = ''
    all_report_files = []
    report_files = []
    
    if 'pipe_steps' in options:
        joined_steps = ', '.join(options['pipe_steps'])
        analysis +=  'Looking For Job Reports for Pipe Steps: %s\n' % joined_steps
        expected_report_files = [(step + '.report') for step in options['pipe_steps']] 

    else:
        analysis += 'Looking for Job Reports.\n'
        expected_report_files = []

    for (path,dirs,files) in os.walk(options['working_directory']):
        for file_name in files:
            if file_name.endswith('.report') and file_name != 'Summary.report':
                all_report_files.append(os.path.join(path, file_name))

    #If expected order for report files, order found files by expected order
    for expected_report_file in expected_report_files:
        for report_file in list(all_report_files):
            if os.path.basename(report_file) == expected_report_file:
               report_files.append(report_file)
               all_report_files.remove(report_file)
               break

    #Add any remaining "unexpected" report files to list.
    report_files += all_report_files

    if report_files:
        analysis += 'Report Files Found:\n'
        for report_file in report_files:
            analysis += ' -- %s\n' % report_file
    else:
        analysis += 'No Report Files Found.\n'
        print analysis
        return analysis

    reports = []

    analysis += '\n'
    for report_file_name in report_files:
        report_file_basename = os.path.basename(report_file_name)
        analysis += report_file_basename + '\n'

        report_file = open(report_file_name, 'r')
        report = report_file.read().strip()
        report_file.close()
        analysis += report + '\n'

        if report_file_basename.endswith('.report'):
            report_name = report_file_basename[:(-len('.report'))]
        else:
            report_name = report_file_basename

        report_type, header, data = read_report(report)
        reports.append((report_name, report_type, header, data))
        analysis += '\n'


    combined_report = combine_report(reports)
    analysis += '\n'
    analysis += 'Combined Report:\n'
    analysis += combined_report + '\n'

    if options['write_report']:
        report_file = os.path.join(options['working_directory'],
                                   JOB_TYPE + '.report')
        write_file(report_file, combined_report)

        if options['write_csv_report']: 
            csv_report = combine_report(reports, write_separator=',')
            report_file = os.path.join(options['working_directory'],
                                       JOB_TYPE + '.report.csv')
            write_file(report_file, csv_report)

    print analysis
    return analysis


def run(options):
    if __name__ != '__main__' and options['is_pipe']:
        out_file = open(options['out_file'], 'w')
        terminal_out, terminal_error = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = out_file, out_file

    print 'Creating Summary Report'

    analyze(options)

    print 'Summary Complete.'                        

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
