#TFLOW Segment: Trim reads using Trimmomatic Read Trimming Utility
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import subprocess
import optparse
import os
import sys
from collections import OrderedDict

if __name__ == "__main__" and __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../'))
    import tflow
    __package__ = "tflow.segments"

from .. import local_settings
if hasattr(local_settings, 'TRIMMOMATIC_LOCATION'):
    TRIMMOMATIC_LOCATION = local_settings.TRIMMOMATIC_LOCATION
else:
    TRIMMOMATIC_LOCATION = ''

if hasattr(local_settings, 'TRIMMOMATIC_EXEC'):
    TRIMMOMATIC_EXEC = local_settings.TRIMMOMATIC_EXEC
else:
    TRIMMOMATIC_EXEC = os.path.join(TRIMMOMATIC_LOCATION, 'trimmomatic-0.32.jar')

from ..util import (print_exit, print_warning, print_error, read_file_list, write_file, 
                    write_file_list, delete_pid_file, count_FASTQ_all, ensure_FASTQ_GZ, 
                    percent_string, print_warning, stop_TFLOW_process)
from .. import util

from .parser_class import OutputParser

JOB_TYPE = 'Trimmomatic'
PROGRAM_URL = 'http://www.usadellab.org/cms/?page=trimmomatic'
SEGMENT_FOR_VERSION = 'trinityrnaseq_r20140717'
COMMAND_LIST = ['java', '-jar', TRIMMOMATIC_EXEC]
COMMAND = ' '.join(COMMAND_LIST)
TEST_COMMAND = '--version'
VERSION = None
VERSION_COMMAND = '--version'
OUT_FILE = 'Trimmomatic.out'
MILESTONES = ['Production Trimmer Started',
              'Writing Trimmed',
              'Trimmomatic Job Done',
              ]
TERMINAL_FLAGS = []
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Exception: ERROR',
                 'Not Found']
DEFAULT_SETTINGS = {'is_paired_reads':True,
                    'max_CPU':'4',
                    #FileNaming
                    'output_directory':'Trimmed_Reads',
                    'left_read_indicator':'R1_',
                    'right_read_indicator':'R2_',
                    'single_reads_list':'single_reads_list',
                    'left_reads_list':'left_read_list',
                    'right_reads_list':'right_read_list',                    
                    'include_unpaired_output':False,
                    #Trim Settings
                    'adapater_trimming':['TruSeq3-PE-2.fa', '1','30','10'],
                    'leading_quality':'30',
                    'trailing_quality':'30',
                    'sliding_window':['25','30'],
                    'minimum_length':'75',
                    #TFLOW_Settings
                    'command':COMMAND,
                    'command_list':COMMAND_LIST,
                    'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    'write_report':True,
                    'write_command':True,
                    'write_pid':True,
                   } 

REQUIRED_SETTINGS = ['is_paired_reads', 'working_directory', 'write_report', 'write_command', 
                     'write_pid']

TRIM_SETTINGS_DICT = OrderedDict()
TRIM_SETTINGS_DICT['adapater_trimming'] = ('ILLUMINACLIP:' + TRIMMOMATIC_LOCATION + '/adapters/')
TRIM_SETTINGS_DICT['leading_quality'] = 'LEADING:'
TRIM_SETTINGS_DICT['trailing_quality'] = 'TRAILING:'
TRIM_SETTINGS_DICT['sliding_window'] = 'SLIDINGWINDOW:'
TRIM_SETTINGS_DICT['minimum_length'] = 'MINLEN:'

class Parser(OutputParser):
    def set_local_defaults(self):
        self.milestones = MILESTONES
        self.terminal_flags = TERMINAL_FLAGS
        self.failure_flags = FAILURE_FLAGS
        self.job_type = JOB_TYPE
        self.tail_length = 30

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
    print '    No Analysis Yet Implemented.'
    return ''

def read(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.read_or_notify()

def stop(options):
    job_pid_file = os.path.join(options['working_directory'],
                                JOB_TYPE + '.auto.pid')
    stop_TFLOW_process(job_pid_file, JOB_TYPE)

def clean(options):
    remove_outfile = (options['mode'] == 'reset')
    out_files = [options['single_reads_list'], options['left_reads_list'], 
                 options['right_reads_list']]

    util.clean_TFLOW_auto_files(options['job_type'], options['project_directory'],
                                options['working_directory'], remove_outfile=remove_outfile, 
                                confirm=options['confirm'], out_files=out_files)

def test(options, silent=False):
    try:
        process = subprocess.Popen(options['command_list'] + [options['test_command']], 
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        process.wait()
        output, error = process.communicate()
        if 'Unable to access jarfile' in output:
            print '%s Cannot Be Found With Shell Command: "%s"' % (JOB_TYPE, options['command'])
            if PROGRAM_URL:
                print 'If Not Installed, %s Can be Downloaded From:\n%s' % (JOB_TYPE, PROGRAM_URL)
            if silent:
                return False

        else:
            if silent:
                return True
            print ' -- %s Found!' % JOB_TYPE

    except subprocess.CalledProcessError as p_exception:
        if int(p_exception.returncode) == 1:
            output =  p_exception.output
            if 'Unable to access jarfile' in output:
                print '%s Cannot Be Found With Shell Command: "%s"' % (JOB_TYPE, options['command'])
                output += ('Error Number: %s\nError Text:\n%s' % (str(p_exception.errno), 
                                                                  p_exception.strerr))
                if silent:
                    return False

            else:
                if silent:
                    return True
                print '%s Found!' % JOB_TYPE

    except OSError as error:
        print '%s Cannot Be Found With Shell Command: "%s"' % (JOB_TYPE, options['command'])
        output = 'Error Number: %s\nError Text:\n%s' % (str(error.errno), error.strerror)

    return output


def run(options):
    if __name__ != '__main__' and options['is_pipe']:
        out_file_stream = open(options['out_file'], 'w')
        terminal_out = sys.stdout
        terminal_error = sys.stderr
        sys.stdout = out_file_stream
        sys.stderr = out_file_stream


    total_count_before = 0
    total_count_after = 0

    if options['is_paired_reads']:
        for option in ['raw_left_reads_list', 'raw_right_reads_list', 'left_read_indicator']:
            if option not in options:
                print_exit('Required Option: %s, for input' % option
                           + ' of paired end reads not found.')

        for reads_list in ['raw_left_reads_list', 'raw_right_reads_list']:
            if not os.path.isfile(options[reads_list]):
                print_exit('Reads List: %s Not Found' % reads_list
                           + ' at Location: $s' % options[reads_list])

        print 'Reading Left Reads File List: %s' % options['raw_left_reads_list'] 
        raw_left_reads = read_file_list(options['raw_left_reads_list'])

        print 'Reading Right Reads File List: %s' % options['raw_right_reads_list'] 
        raw_right_reads = read_file_list(options['raw_right_reads_list'])

        if len(raw_left_reads) != len(raw_right_reads):
            print_exit('Number of Left Reads: %i' % len(raw_left_reads)
                       + 'Does Not Equal Number of Right Reads: %i' % len(raw_right_reads))

        raw_reads = raw_left_reads + raw_right_reads

    else:
        if not 'raw_single_reads_list' in options:
            print_exit('Required Option: raw_single_reads_list, for input'
                       + ' of unpaired reads not found.')

        if not os.path.isfile(options['raw_single_reads_list']):
            print options['raw_single_reads_list']
            print_exit('Reads List: %s Not Found' % 'raw_single_reads_list'
                       + ' at Location: %s' % options['raw_single_reads_list'])

        print 'Reading Unpaired Reads File List: %s' % options['raw_single_reads_list'] 
        raw_reads = read_file_list(options['raw_single_reads_list'])

    for read in raw_reads:
        full_read = os.path.join(options['working_directory'], read)
        ensure_FASTQ_GZ(full_read)

    print ''
    print 'Production Trimmer Started, going to trim files:'
    if options['is_paired_reads']:
        print '  Left Reads:'
        for read in raw_left_reads:
            print('  -- ' + read)
        print ''
        print '  Right Reads:'
        for read in raw_right_reads:
            print('  -- ' + read)
        print ''
        trim_reads = raw_left_reads

    else:
        print '  Unpaired Reads:'
        for read in raw_reads:
            print('  -- ' + read)
        print ''
        trim_reads = raw_reads

    if options['write_command']:
            command_file_name = os.path.join(options['project_directory'], 
                                             options['job_type'] + '.auto.sh')
            command_file = open(command_file_name, 'w')
            command_file.write('#!/bin/sh\n')

    if 'output_directory' in options:
        full_output_directory = os.path.join(options['working_directory'],
                                             options['output_directory'])
        if not os.path.isdir(full_output_directory):
            print 'Preparing Output Directory: %s' % full_output_directory
            os.makedirs(full_output_directory)

    out_files = OrderedDict()
    if options['is_paired_reads']:
        out_files['left_paired'] = []
        out_files['right_paired'] = []
        out_files['left_unpaired'] = []
        out_files['right_unpaired'] = []
    else:
        out_files['single_reads'] = []

    for read in trim_reads:
        full_read = os.path.join(options['working_directory'], read)
        if not os.path.isfile(full_read):
            print_exit('Input Read File for Trimming: %s Cannot Be Found.' % full_read)

        starting_reads = count_FASTQ_all(full_read)
        total_count_before += starting_reads

        print ''
        print ' ----- Trimming File: %s  -----' % read
        print ''
        print '  Starting Reads:', starting_reads

        threads_list = ['-threads', options['max_CPU']]
        threads = ' '.join(threads_list)

        #Prepare Command for Paired Reads
        if options['is_paired_reads']:
            mode = 'PE'
            base_name = os.path.basename(read)
            base_name = base_name.rstrip('gz').rstrip('.').rstrip('fastq').rstrip('.')

            out_base_name = ''.join(base_name.split(options['left_read_indicator']))
            out_read = out_base_name + '.fq'

            if out_base_name == base_name:
                out_base_name = base_name + '-Trimmed'
                out_read = out_base_name + '.fq'
                print ('Identifier: "%s" for Left Read ' % options['left_read_indicator']
                       + 'Not Found. Using Default Naming: %s' % out_read)

            if 'output_directory' in options:
                out_read = os.path.join(options['output_directory'], out_read)
                out_base_name = os.path.join(options['output_directory'], out_base_name)

            full_out_read = os.path.join(options['working_directory'], out_read)

            base_out_list = ['-baseout', full_out_read]
            base_out = ' '.join(base_out_list)
            base_in_list = ['-basein', full_read]
            base_in = ' '.join(base_in_list)
            trim_log_name =  base_name + '.trimlog'
            if 'output_directory' in options:
                trim_log_name = os.path.join(options['output_directory'], trim_log_name)
            trim_log_list = ['-trimlog', trim_log_name]
            trim_log = ' '.join(trim_log_list)


            expected_out_files = OrderedDict()
            expected_out_files['left_paired'] = out_base_name + '_1P.fq'
            expected_out_files['right_paired'] = out_base_name + '_2P.fq'
            expected_out_files['left_unpaired'] = out_base_name + '_1U.fq'
            expected_out_files['right_unpaired'] = out_base_name + '_2U.fq'

            expected_full_out_files = OrderedDict()
            for out_file in expected_out_files:
                expected_full_out_files[out_file] = os.path.join(options['working_directory'],
                                                                 expected_out_files[out_file])

            print ''
            print '  Output File Basename: %s' % out_base_name
            print '  Expected Output Files:'
            for file_name in expected_out_files.values():
                print '   ', file_name

            action_list = []
            print ''
            print '  Trim Settings:'
            for trim_setting in TRIM_SETTINGS_DICT:
                if trim_setting in options:
                    setting = options[trim_setting]
                    if isinstance(setting, list):
                        setting = ':'.join(options[trim_setting])

                    formatted_setting = (TRIM_SETTINGS_DICT[trim_setting] + setting)
                    action_list.append(formatted_setting)
                    print '   ', formatted_setting.strip()

            command_segments = ([options['command'], mode, threads,
                                 base_out, base_in, trim_log] + action_list)

            command_list = (list(options['command_list']) + [mode] + threads_list + base_out_list 
                            + base_in_list + trim_log_list + action_list)

        #Prepare Command for Unpaired Reads
        else:
            mode = 'SE'
            base_name = os.path.basename(read)
            base_name = base_name.rstrip('gz').rstrip('.').rstrip('fastq').rstrip('.')

            out_base_name = base_name + '-Trimmed'
            out_read = out_base_name + '.fq'

            if 'output_directory' in options:
                out_read = os.path.join(options['output_directory'], out_read)
                out_base_name = os.path.join(options['output_directory'], out_base_name)

            full_out_read = os.path.join(options['working_directory'], out_read)

            trim_log_name =  base_name + '.trimlog'
            if 'output_directory' in options:
                trim_log_name = os.path.join(options['output_directory'], trim_log_name)
            trim_log_list = ['-trimlog', trim_log_name]
            trim_log = ' '.join(trim_log_list)


            expected_out_files = OrderedDict()
            expected_out_files['single_reads'] = out_read

            expected_full_out_files = OrderedDict()
            for out_file in expected_out_files:
                expected_full_out_files[out_file] = os.path.join(options['working_directory'],
                                                                 expected_out_files[out_file])

            print ''
            print '  Output File Basename: %s' % out_base_name
            print '  Expected Output Files:'
            for file_name in expected_out_files.values():
                print '   ', file_name

            action_list = []
            print ''
            print '  Trim Settings:'
            for trim_setting in TRIM_SETTINGS_DICT:
                if trim_setting in options:
                    setting = options[trim_setting]
                    if isinstance(setting, list):
                        setting = ':'.join(options[trim_setting])

                    formatted_setting = (TRIM_SETTINGS_DICT[trim_setting] + setting)
                    action_list.append(formatted_setting)
                    print '   ', formatted_setting.strip()

            command_segments = ([options['command'], mode, threads, trim_log,
                                 full_read, full_out_read] + action_list)

            command_list = (list(options['command_list']) + [mode] + threads_list + trim_log_list
                            + [full_read] + [full_out_read] + action_list)


        #Run Trimmomatic Command
        print ''  
        print '  Running Command with Segments:'
        for segment in command_segments:
            print '    ' + segment

        command = ' '.join(command_list)
        if options['write_command']:
            command_file.write(command + '\n\n')

        print ''
        print command
        print ''
        sys.stdout.flush()

        try:
            process = subprocess.Popen(command_list, stdout=sys.stdout, stderr=sys.stderr,
                                       cwd=options['working_directory'])
            if options['write_pid']:
                pid_file_name = os.path.join(options['working_directory'], 
                                             options['job_type'] + '.auto.pid')
                write_file(pid_file_name, str(process.pid))
            process.wait()

            if options['write_pid']:
                delete_pid_file(pid_file_name)

        except KeyboardInterrupt:
            if __name__ != '__main__' and options['is_pipe']:
                sys.stdout, sys.stderr = terminal_out, terminal_error
                out_file_stream.close()
            print ''
            print 'Killing Trimmomatic Process'
            process.kill()
            raise

        sys.stdout.flush()

        print ''
        print '  Output Files:'
        if not any(os.path.isfile(x) for x in expected_full_out_files.values()):
            print ' ', 
            print_warning('No Output Files Found for File: %s ' % read )
            print ''
            print ''
            continue

        else:
            for out_file in expected_full_out_files:
                if os.path.isfile(expected_full_out_files[out_file]):
                    print '    Found:', expected_out_files[out_file]
                    print out_files
                    print expected_out_files
                    out_files[out_file].append(expected_out_files[out_file])

        if options['is_paired_reads']:
            final_reads = count_FASTQ_all(os.path.join(options['working_directory'],
                                                       expected_out_files['left_paired']))
        else:
            final_reads = count_FASTQ_all(os.path.join(options['working_directory'],
                                                       expected_out_files['single_reads']))

        total_count_after += final_reads
        print ''
        print 'Finished With File:', read
        print 'Starting Reads:', starting_reads
        print 'Final Reads:', final_reads
        print percent_string(final_reads, starting_reads), 'Remaining.'
        print ''

    command_file.close()

    print 'Initial Read Count:', total_count_before
    print 'Final Read Count:  ', total_count_after
    print percent_string(total_count_after, total_count_before), 'of Original'

    print ''
    print 'Writing Final Output Files:'
    if options['is_paired_reads']:
        left_reads = list(out_files['left_paired'])
        right_reads = list(out_files['right_paired'])
        if options['include_unpaired_output']:
            left_reads += list(out_files['left_unpaired'])
            right_reads += list(out_files['right_unpaired'])
        if not left_reads:
            print_exit('No Ouput Left Reads Were Found!')
        elif not right_reads:
            print_exit('No Ouput Right Reads Were Found!')
        print 'Writing Trimmed Left Read Files to List: %s' % options['left_reads_list']
        for read in left_reads:
            print ' ', read
        print ''
        print 'Writing Trimmed Right Read Files to List: %s' % options['right_reads_list']
        for read in right_reads:
            print ' ', read

        write_file_list(os.path.join(options['working_directory'], options['left_reads_list']),
                        left_reads)

        write_file_list(os.path.join(options['working_directory'], options['right_reads_list']),
                        right_reads)

    else:
        single_reads = out_files['single_reads']
        if not single_reads:
            print_exit('No Ouput Trimmed Reads Were Found!')
        print 'Writing Trimmed Read Files to List: %s' % options['single_reads_list']
        for read in single_reads:
            print ' ', read

        write_file_list(os.path.join(options['working_directory'], options['single_reads_list']),
                        single_reads)

    print 'Trimmomatic Job Done.'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout = terminal_out
        sys.stderr = terminal_error
        out_file_stream.close()


