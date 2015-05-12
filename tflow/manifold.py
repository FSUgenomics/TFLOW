#!/usr/bin/env python2.7
#TFLOW Component: General Job Tracker
#Usage: "manifold.py RUN_MODE [--options]"
#For Advanced Usage: "manifold.py -h"
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import sys
import os
import argparse
from copy import deepcopy

if __name__ == "__main__" and __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
    import tflow
    __package__ = "tflow"

from .util import (get_file_settings, print_warning, print_except, print_multi, print_exit,
                   write_file, print_settings, write_settings, write_date_time, lowercase, 
                   flexible_boolean_string, BOOL, FLEXIBLE_BOOL, ACTION_NAMES, DEFAULT_SETTINGS)

MODES = ['track', 'analyze', 'run', 'read', 'test', 'stop', 'print_settings']
NULL_OUT_FILES = ['N/A', 'n/a', 'NA', 'na', 'None', 'none', None, '']
READ_TYPES = ['fq', 'fa']
JOB_TYPES = []

def parse_args():
    parser = argparse.ArgumentParser(prog='manifold.py', 
                                     description='Run, Track, and Analyze Assembly Jobs.')
    tflow_args = parser.add_argument_group('General TFLOW Args', 'Arguments for TFLOW Runs.')
    tflow_args.add_argument('mode', action='store', default='track', type=lowercase, 
                            help='Select Conductor Mode', choices=MODES, metavar='MODE')
    #parser.add_argument('-m', '--mode', action='store', default='track', type=lowercase, 
    #			help='Select Conductor Mode', choices=MODES, metavar='MODE')
    tflow_args.add_argument('-t', '--job_type', action='store',  
                            help='Job Type to Conduct', metavar='JOB_TYPE')
    tflow_args.add_argument('-o', '--out_file', action='store', 
                            help='Job Output File to Evaluate', metavar='OUT_FILE')
    tflow_args.add_argument('-v', '--verbose', action='store_true', default=None,
                            help='Verbose Line Output', dest='verbose')
    tflow_args.add_argument('--overwrite', action='store', default=None, 
                            type=flexible_boolean_string, help='Overwrite Previous Outputs', 
                            choices=BOOL, metavar='BOOL')
    tflow_args.add_argument('--is_paired_reads', action='store', default=None, 
                            type=flexible_boolean_string,
                            help='If Working With Reads, Whether Are Paired-End Reads',
                            choices=BOOL, metavar='BOOL')
    tflow_args.add_argument('--read_type', action='store', default=None, 
                            help='If Working With Reads, Type of Reads. (FASTQ or FASTA)',
                            choices=READ_TYPES, metavar='TYPE')
    tflow_args.add_argument('--label', action='store', default=None, 
                            help='Tissue-Type Label for Automated Labeling', metavar='LABEL')
    tflow_args.add_argument('--max_CPU', action='store', default=None, 
                            help='Maximum Number of CPU\'s to Use for Run.', metavar='#CPU')

    testing_args = parser.add_argument_group('Test Mode Args', 
                                              'Arguments for Test Mode')
    testing_args.add_argument('--print_test_output', action='store', default=None, 
                            type=flexible_boolean_string, help='Print Output of Test Commands', 
                            choices=BOOL, metavar='BOOL')

    trimming_args = parser.add_argument_group('Read Trimming Args', 
                                              'Arguments for "Make Reads List and '
                                              + 'Read Trimming Before Assembly')
    trimming_args.add_argument('--raw_reads', action='store', nargs='*', default=None,
                               help='Reads Files For Trimming')
    trimming_args.add_argument('--raw_left_reads', action='store', nargs='*', default=None,
                               help='Forward Paired-End Reads Files For Trimming')
    trimming_args.add_argument('--raw_right_reads', action='store', nargs='*', default=None,
                               help='Reverse Paired-End Reads Files For Trimming')
    trimming_args.add_argument('--raw_reads_list', action='store', default=None,
                               help='List of Reads Files For Trimming')
    trimming_args.add_argument('--raw_left_reads_list', action='store', default=None,
                               help='List of Forward Paired-End Reads Files For Trimming')
    trimming_args.add_argument('--raw_right_reads_list', action='store', default=None,
                               help='List of Reverse Paired-End Reads Files For Trimming')
    trimming_args.add_argument('--include_unpaired_output', action='store', default=None,
                               help='Include Unpaired Reads in Final Read Lists')
    trimming_args.add_argument('--left_read_indicator', action='store', default=None,
                               help='Indicator of left reads, for automated parsing when paired')
    trimming_args.add_argument('--right_read_indicator', action='store', default=None,
                               help='Indicator of right reads, for automated parsing when paired')

    assembly_args = parser.add_argument_group('Assembly Args', 
                                              'Arguments for Assembly of Reads/Sequences')
    assembly_args.add_argument('--reads', action='store', nargs='*', default=None,
                               help='Reads Files for Assembly')
    assembly_args.add_argument('--left_reads', action='store', nargs='*', default=None,
                               help='Forward Paired-End Reads Files for Assembly')
    assembly_args.add_argument('--right_reads', action='store', nargs='*', default=None,
                               help='Reverse Paired-End Reads Files for Assembly')
    assembly_args.add_argument('--reads_list', action='store', default=None,
                               help='List of Reads Files For Assembly')
    assembly_args.add_argument('--left_reads_list', action='store', default=None,
                               help='List of Forward Paired-End Reads Files For Assembly')
    assembly_args.add_argument('--right_reads_list', action='store', default=None,
                               help='List of Reverse Paired-End Reads Files For Assembly')
    assembly_args.add_argument('--relative_input_file', action='store', default=None,
                               help='Input File for CAP3 Assembly, Relative to Project Directory')
    assembly_args.add_argument('--absolute_input_file', action='store', default=None,
                               help='Absolute Path to Input File for CAP3 Assembly')
    assembly_args.add_argument('--relative_input_files', action='store', nargs = '*',
                               default=None, help=('Input Files for CAP3 Assembly, Relative to'
                                                   + ' Project Directory'))
    assembly_args.add_argument('--absolute_input_files', action='store', nargs='*', default=None,
                               help='Absolute Path to Input Files for CAP3 Assembly')

    analysis_args = parser.add_argument_group('Analysis Args', 
                                              'Arguments for Analysis of Sequence Files')

    analysis_args.add_argument('--rel_input_analysis_file', action='store', default=None,
                               help='Input File for Analysis, Relative to Project Directory')
    analysis_args.add_argument('--absolute_input_analysis_file', action='store', default=None,
                               help='Absolute Path to Input File for Analysis')
    analysis_args.add_argument('--BUSCO_type', action='store', default=None,
                               help='Organism Type for BUSCO Analysis.')
    analysis_args.add_argument('--copy_input_file', action='store', default=None, 
                               type=flexible_boolean_string, help='Copy Source Input File Instead of Using Inplace.',
                               choices=BOOL, metavar='BOOL')

    #parser.set_defaults()

    return (vars(parser.parse_args()))


def get_settings():
    args = parse_args()
    settings = get_file_settings()
    for arg_name in args:
        if args[arg_name] not in ['', None]:
            settings[arg_name]=args[arg_name]
    
    for setting in DEFAULT_SETTINGS:
        if setting not in settings:
            settings[setting] = DEFAULT_SETTINGS[setting]

    required = {'job_type':'Job Type Not Specified',
                }
                
    for required_arg in required:
        if required_arg not in settings:
            print required[required_arg]
            sys.exit(1)           

    segments_module = __import__('tflow.segments', fromlist=[settings['job_type']])
    pipes_module = __import__('tflow.pipes', fromlist=[settings['job_type']])
    if (not hasattr(segments_module, settings['job_type'])
        and not hasattr(pipes_module, settings['job_type'])):
        print_exit('Job Type: %s Not Found.' % settings['job_type'], 1)
    
    return settings


def flow(options, check_done=False):
    job_type = options['job_type']
    segments_module = __import__('tflow.segments', fromlist=[job_type])
    module = getattr(segments_module, job_type)

    job_options = deepcopy(options)

    if hasattr(module, 'DEFAULT_SETTINGS'):
        default_settings = module.DEFAULT_SETTINGS
        for setting in default_settings:
            if setting not in job_options:
                job_options[setting] = default_settings[setting]

    if 'out_file' in options:
        out_file = options['out_file']
    else:
        if hasattr(module, 'OUT_FILE'):
            out_file = module.OUT_FILE
        else:
            print_except('For Job Type %s, No Output File Given.' % job_type)

    if 'working_directory' in options and options['working_directory'] not in [None, 'None', 'none', 'N/A', 'n/a']:
        working_directory = options['working_directory']
        if out_file:
            out_file = os.path.join(working_directory, out_file)
        if not os.path.isdir(working_directory):
            if options['mode'] == 'run':
                full_working_directory = os.path.join(options['project_directory'],
                                                      options['working_directory'])
                print 'Making New Working Directory for Run: %s' % options['working_directory']
                os.makedirs(full_working_directory)
            elif options['mode'] not in ['test']:
                print_warning('Working Directory %s Not Found.' % working_directory)

    else:
        working_directory = ''

    if (out_file and os.path.basename(out_file) not in NULL_OUT_FILES 
        and not os.path.isfile(out_file) and job_options['mode'] not in ['run', 'test']):
        print_exit(['Job Output File %s Not Found...' % out_file], 1)

    job_options['out_file'] = out_file
    terminal_output = (sys.stdout, sys.stderr)

    if check_done:
        if not hasattr(module, 'check_done'):
            print_except('Job Type %s Has No check_done Method.' % job_type)

        try:
            return module.check_done(job_options)

        except KeyboardInterrupt:
            (sys.stdout, sys.stderr) = terminal_output
            print_exit('')

    elif options['mode'] == 'run':
        print_multi('', 'Running %s Job...' % job_type, '')
        if not hasattr(module, 'run'):
            print_except('Job Type %s Has No Run Method.' % job_type)

        try:
            if options['write_settings']:
                settings_file_name = os.path.join(options['working_directory'], 
                                                  options['job_type'] + '.auto.settings')
                write_settings(job_options, settings_file_name)

            if options['write_times']:
                time_file = os.path.join(options['working_directory'], 
                                         options['job_type'] + '.auto.timing')
                start_time = write_date_time(time_file)

            module.run(job_options)

            if options['write_times']:
                write_date_time(time_file, start=start_time)

        except KeyboardInterrupt:
            (sys.stdout, sys.stderr) = terminal_output
            print_exit(['', 'Running Stopped.'], 2)

    elif options['mode'] == 'track':
        print_multi('', 'Tracking %s Job...' % job_type, '')
        if not hasattr(module, 'track'):
            print_except('Job Type %s Has No Track Method.' % job_type)

        try:
            module.track(job_options)

        except KeyboardInterrupt:
            (sys.stdout, sys.stderr) = terminal_output
            print_exit(['', 'Tracking Stopped.'], 2)

    elif options['mode'] == 'analyze':
        print_multi('Analyzing %s Job...' % job_type, '')
        if not hasattr(module, 'analyze'):
            print_except('Job Type %s Has No Analyze Method.' % job_type)

        try:
            analysis = module.analyze(job_options)

        except KeyboardInterrupt:
            (sys.stdout, sys.stderr) = terminal_output
            print_exit(['', 'Analysis Stopped.'], 2)

        if analysis and options['write_analysis']:
            analysis_file_name = os.path.join(options['working_directory'],
                                              options['job_type'] + '.auto.analysis')
            write_file(analysis_file_name, analysis)

    elif options['mode'] == 'read':
        print_multi('Reading %s Job...' % job_type, '')
        if not hasattr(module, 'read'):
            print_except('Job Type %s Has No Read Method.' % job_type)

        try:
            analysis = module.read(job_options)

        except KeyboardInterrupt:
            (sys.stdout, sys.stderr) = terminal_output
            print_exit(['', 'Reading Stopped.'], 2)

    elif options['mode'] == 'test':
        print_multi('Testing %s Job...' % job_type, '')
        if not hasattr(module, 'test'):
            print_except('Job Type %s Has No Test Method.' % job_type)

        try:
            output = module.test(job_options)
            if options['print_test_output']:
                print output

        except KeyboardInterrupt:
            (sys.stdout, sys.stderr) = terminal_output
            print_exit(['', 'Testing Stopped.'], 2)


    print ''


def segment(options):
    #Add Task-Specific Options From Options File to Segment Settings
    job_type = options['job_type']   
    if job_type in options:
        job_dict = options[job_type]
        for setting in job_dict:
            if setting in options:               
                print_warning('Option:  "%s"  with' % setting
                              + ' value: "%s"  Being Overridden' % options[setting]
                              + ' for job  "%s" ' % job_type
                              + ' by Job-Specific Options File'
                              + ' Setting:  "%s"' % job_dict[setting])

            options[setting] = job_dict[setting]

    #Remove Unused Step Settings from Segment Settings
    for item in options.keys():
        if isinstance(options[item], dict):
            del(options[item])

    if 'working_directory' not in options:
        options['working_directory'] = options['project_directory']

    step_done = flow(options, check_done=True)

    if options['mode'] == 'run':
        if step_done and not options['overwrite']:
            print 'Running %s Job.\n\n    %s Job Already Complete.' % (job_type, job_type)
        else:
            flow(options)
            print '    %s Job Complete.' % job_type

    elif options['mode'] == 'track':
        if not step_done:
            flow(options)
        else:
            print'Tracking %s Job.\n' % job_type
        print '    %s Job Complete.' % job_type

    elif options['mode'] == 'analyze':
        if not step_done:
            print ('Analyzing %s Job.\n\n%s Job ' % (job_type, job_type)
                       + 'Not Yet Complete, Cannot Analyze.')
        else:
            flow(options)
            print '    %s Analysis Complete.' % job_type

    elif options['mode'] in ['test', 'read', 'stop']:
        flow(options)

    else:
        print_except('Conductor Has Unrecognized Mode Type %s' % options['mode'])
 
    print ''


def manifold(options):
    #Get Segment Pipe Module Object
    segments_module = __import__('tflow.pipes', fromlist=[options['job_type']])
    module = getattr(segments_module, options['job_type'])
  
    #Get Steps from Pipe Module
    if hasattr(module, 'steps'):
        pipe_steps = getattr(module, 'steps')
    else:
        print_except('Orchestrator Module: %s  Does not Have Steps Parameter.' % pipe)

    #Print Steps
    print_multi('%s Steps: %s in Pipe: %s' % (ACTION_NAMES[options['mode']], 
                                              ', '.join(pipe_steps.keys()), 
                                              options['job_type']), '')

    #Test for Existence of Each Step
    for step in pipe_steps:
        step_module = __import__('tflow.segments', fromlist=[step])
        if not hasattr(step_module, step):
            print_exit('Job Type: %s Not Found.' % step, 1)

    #Add List of steps to options.
        options['pipe_steps'] = list(pipe_steps)

    #If Running, Write Timing and Settings
    if options['mode'] == 'run':
        if options['write_settings']:
            pipe_settings_file_name = os.path.join(options['project_directory'], 
                                                   options['job_type'] + '.auto.settings')
            write_settings(options, pipe_settings_file_name)

        if options['write_times']:
            pipe_time_file = os.path.join(options['project_directory'], 
                                     options['job_type'] + '.auto.timing')
            pipe_start_time = write_date_time(pipe_time_file)

    #Perform Each Step
    for step in pipe_steps:
        step_run_options = deepcopy(options)
        pipe_step_options = pipe_steps[step]

        #Remove Step Settings from Segment Settings
        for item in step_run_options.keys():
            if isinstance(step_run_options[item], dict):
                del(step_run_options[item])

        #Add Task-Specific Options From Options File to Segment Settings
        if step in options:
            step_dict = options[step]
            for setting in step_dict:
                if setting in step_run_options:
                    print_warning('Option  "%s" ' % setting
                                  + ' for step  "%s" ' % step
                                  + ' value "%s"' % str(step_run_options[setting])
                                  + ' is being overwritten by step-specific options file'
                                  + ' value:  "%s"' % step_run_options[setting])
                step_run_options[setting] = step_dict[setting]


        #Add Pipe-Specific Step Settings to Segment Settings
        for setting in pipe_step_options:
            #If pipe-specific option already given, print warning and override.
            if setting in step_run_options:
                print_warning('Option  "%s" ' % setting
                              + ' with Value:  "%s" ' % step_run_options[setting]
                              + ' is Being Overridden for Pipe Step  "%s" ' % step
                              + ' by Pipe Setting Value:  "%s" ' % pipe_step_options[setting])
            step_run_options[setting] = pipe_step_options[setting]

        if 'working_directory' in pipe_step_options:
            full_working_directory = os.path.join(options['project_directory'], 
                                                  pipe_step_options['working_directory']) 
            step_run_options['working_directory'] = full_working_directory 

        else:
           step_run_options['working_directory'] = options['project_directory'] 

        step_run_options['job_type'] = step

        step_done = flow(step_run_options, check_done=True)

        if options['mode'] == 'run':
            if step_done and not options['overwrite']:
                print 'Running %s Job.\n\n    %s Job Already Complete.' % (step, step)
            else:
                flow(step_run_options)
                print '    %s Job Complete.' % step

        elif options['mode'] == 'track':
            if not step_done:
                flow(step_run_options)
            else:
                print'Tracking %s Job.\n' % step
            print '    %s Job Complete.' % step

        elif options['mode'] == 'analyze':
            if not step_done:
                print ('Analyzing %s Job.\n\n%s Job ' % (step, step)
                       + 'Not Yet Complete, Cannot Analyze.')
            else:
                flow(step_run_options)
                print '%s Analysis Complete.' % step

        elif options['mode'] in ['test', 'read', 'stop']:
            flow(step_run_options)
 
        else:
            print_except('Conductor Has Unrecognized Mode Type %s' % options['mode'])

        print ''
        sys.stdout.flush()

    #Write End Time of Pipe
    if options['mode'] == 'run' and options['write_times']:
        write_date_time(pipe_time_file, start=pipe_start_time)




if __name__ ==  '__main__':

    options = get_settings()    

    if 'project_directory' not in options:
        options['project_directory'] = os.getcwd()

    print ''

    if options['mode'] == 'print_settings':
        print_settings(options, 'TFLOW Conductor Options:')
        print ''
        sys.exit(0)

    if any(x in options['job_type'] for x in ['_pipe', '_Pipe']):
        options['is_pipe'] = True
        manifold(options)
    else:
        segment(options)

    print 'All Jobs Complete.'
    print ''
