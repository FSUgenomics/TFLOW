#TFLOW Segment: Find Gene Annotations for Query Sequence File useing Reference Nucleotide or 
#               Protein Sequence File
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 05/22/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import os.path
import sys
import subprocess
import shutil

REFERENCE_TYPES = {'Protein':'prot', 'protein':'prot', 'prot':'prot', 
                   'Nucleotide':'nucl', 'nucleotide':'nucl', 'nucl':'nucl'}

if __name__ == "__main__" or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../'))
    import tflow.segments
    __package__ = "tflow.segments"

from .parser_class import OutputParser
from ..util import print_exit, print_return
from .. import util
from .. import local_settings
from ..annotation import (Annotation, Annotation_Record, Annotation_Database, 
                          Annotation_Map, Name_Map)

if hasattr(local_settings, 'BLAST_LOCATION'):
    BLAST_LOCATION = local_settings.BLAST_LOCATION
else:
    BLAST_LOCATION = ''

if hasattr(local_settings, 'BLAST_EXEC'):
    BLAST_EXEC = local_settings.BLAST_EXEC
else:
    BLAST_EXEC = os.path.join(BLAST_LOCATION, 'blastx')

if hasattr(local_settings, 'MAKE_BLAST_DB_LOCATION'):
    MAKE_BLAST_DB_LOCATION = local_settings.MAKE_BLAST_DB_LOCATION
else:
    MAKE_BLAST_DB_LOCATION = ''

if hasattr(local_settings, 'MAKE_BLAST_DB_EXEC'):
    MAKE_BLAST_DB_EXEC = local_settings.MAKE_BLAST_DB_EXEC
else:
    MAKE_BLAST_DB_EXEC = os.path.join(MAKE_BLAST_DB_LOCATION, 'makeblastdb')


JOB_TYPE = 'Prototype_Find_Annotations'
PROGRAM_URL = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download'
SEGMENT_FOR_VERSION = '2.2.29'
BLAST_COMMAND = BLAST_EXEC
BLAST_COMMAND_LIST = [BLAST_COMMAND]
BLAST_DB_COMMAND = MAKE_BLAST_DB_EXEC
BLAST_DB_COMMAND_LIST = [BLAST_DB_COMMAND]
TEST_COMMAND = '-h'
OUT_FILE = 'Prototype_Find_Annotations.out'
MILESTONES = ['Annotation Complete']
TERMINAL_FLAGS = []
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Exception: ERROR',
                 'Not Found',
                 ]

MATCH_PREFIX = 'Matches'
ANNOTATION_PREFIX = 'Annotations'

DEFAULT_SETTINGS = {'copy_input_file':False,
                    'max_CPU':'4',
                    'blast_evalue':'1e-5',
                    'evalue_cutoffs':['1e-10', '1e-20', '1e-40'],
                    'blast_result_file':'blast.out',
                    'db_title':'BLAST_DB',
                    'max_matches':'5',
                    'write_all_matches':True,
                    'write_best_matches':True,
                    'verbose_tracking':True,
                    #'reference_type':'nucl',
                    #TFLOW Settings
                    'blast_command':BLAST_COMMAND,
                    'blast_command_list':BLAST_COMMAND_LIST,
                    'blast_db_command':BLAST_DB_COMMAND,
                    'blast_db_command_list':BLAST_DB_COMMAND_LIST,
                    'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    'write_report':True,
                    'write_command':True,
                    'write_pid':True,
                    }

REQUIRED_SETTINGS = ['blast_command_list', 'blast_db_command_list', 'working_directory', 
                     'copy_input_file', 'blast_evalue', 'max_CPU', 'blast_result_file', 
                     'evalue_cutoffs', 'write_command', 'write_report', 'write_pid', 
                     'reference_type', 'max_matches']

REQUIRED_ANALYSIS_SETTINGS = ['blast_result_file', 'evalue_cutoffs', 'write_report', 
                              'write_all_matches', 'write_best_matches']

OPTIONAL_TRACKING_SETTINGS = ['working_directory', 'blast_result_file']


class Parser(OutputParser):
    def set_local_defaults(self):
        self.milestones = MILESTONES
        self.terminal_flags = TERMINAL_FLAGS
        self.failure_flags = FAILURE_FLAGS
        self.job_type = JOB_TYPE

    def check_queries_processed(self, blast_file_name):
        count = 0
        if self.shell_tracking:
            try:
                string_count = subprocess.check_output(['grep', '-c', '# Query:', 
                                                        blast_file_name])
                count = int(string_count.strip())
            except subprocess.CalledProcessError:
                util.print_warning('Problem with Shell Query Checking? Continuing...')
        else:
            with open(blast_file_name, 'r') as blast_file:
                for line in blast_file:
                    if line.startswith('# Query:'):
                        count += 1
        return count

    def annotation_track(self, options, loud=False):
        from time import sleep
        if 'verbose_tracking' in options and options['verbose_tracking']:
            verbose_tracking = True
            #Ensure Required Settings in Options
            for optional_option in OPTIONAL_TRACKING_SETTINGS:
                if optional_option not in options:
                    print 'Optional Option: %s for %s tracking not given.' % (optional_option, 
                                                                              JOB_TYPE)
                    print 'Defaulting to non-verbose tracking.'
                    verbose_tracking = False
        else:
            verbose_tracking = False
        
        if verbose_tracking:
            print 'Attempting to Initiate Verbose Annotation Tracking:'
            print ''

        #Ensure Working Directory Exists
        if verbose_tracking and not os.path.isdir(options['working_directory']):
            print 'Working Directory: %s Not Found.' % options['working_directory']
            print 'Defaulting to non-verbose tracking.'
            verbose_tracking = False

        #Ensure A Type of Query Sequence File is Given
        if verbose_tracking and not any(x in options for x in ['absolute_input_analysis_file',
                                                               'rel_input_analysis_file',
                                                               'input_analysis_file',
                                                               'result_name_file']):
            print ('Either input_analysis_file, absolute_input_analysis_file,'
                   + ' rel_input_analysis_file, or result_name_file paramater required.')
            print 'Defaulting to non-verbose tracking.'
            verbose_tracking = False

        #Assign Correct Input File Name
        if verbose_tracking:
            if 'input_analysis_file' in options:
                if os.path.isabs(options['input_analysis_file']):
                    full_input_file = options['input_analysis_file']
                else:
                    full_input_file = os.path.join(options['working_directory'], 
                                                   options['input_analysis_file'])
                                                   
            elif 'absolute_input_analysis_file' in options:
                full_input_file = options['absolute_input_analysis_file']

            elif 'rel_input_analysis_file' in options:
                full_input_file = os.path.join(options['project_directory'], 
                                               options['rel_input_analysis_file'])

            elif 'result_name_file' in options:
                full_result_name_file = os.path.join(options['project_directory'], 
                                                     options['result_name_file'])
                if not os.path.isfile(full_result_name_file):
                    print('Provided File: %s Containing' % full_result_name_file  
                          + ' Result Sequence File Name Not Found.')
                    print 'Defaulting to non-verbose tracking.'
                    verbose_tracking = False

                if verbose_tracking:
                    rf = open(full_result_name_file, 'r')
                    full_input_file = rf.read().strip()
                    rf.close()

                    if not os.path.isfile(full_input_file):
                        print('Cannot Find Read Result Sequence File: %s' % full_input_file)
                        print 'Defaulting to non-verbose tracking.'
                        verbose_tracking = False

        if verbose_tracking:
            input_file = os.path.basename(full_input_file)
            if not any(x(full_input_file) for x in [util.is_FASTA, util.is_FASTA_GZ]):
                print 'File: %s is not a FASTA File.' % full_input_file
                print 'Defaulting to non-verbose tracking.'
                verbose_tracking = False

        if verbose_tracking:
            #Set BLAST Output File
            blast_file = options['blast_result_file']
            full_blast_file = os.path.join(options['working_directory'], blast_file)

            #Check that Blast File Exists
            if not os.path.isfile(full_blast_file):
                print 'Input Reference Sequence File: %s Not Found.' % full_blast_file
                print 'Defaulting to non-verbose tracking.'
                verbose_tracking = False

        if verbose_tracking:
            num_query_sequences = util.count_FASTA_all(full_input_file)
            print 'Verbose Tracking Successfully Initiated.'
            print 'Query Sequences:', num_query_sequences
        else:
            num_query_sequences = 0

        num_queries_processed = 0
        new_queries_processed = 0
        while self.running:
            if verbose_tracking:
                new_queries_processed = self.check_queries_processed(full_blast_file)
            if self.check_updated() or num_queries_processed != new_queries_processed:
                #print 'Updated:', self.check_updated()
                #print 'Num_Queries_New:',  num_queries_processed != new_queries_processed
                self.running = self.check(loud)
                if verbose_tracking:
                    print new_queries_processed, '/', num_query_sequences, 
                    print '(%s)' % util.percent_string(new_queries_processed, 
                                                       num_query_sequences),
                    print 'Query Sequences Processed.'
                    num_queries_processed = new_queries_processed
            if self.running:
                sleep(self.sleep_time)


def check_done(options):
    parser = Parser()
    parser.out_file = options['out_file']
    failure_exit = (options['mode'] in ['run', 'track'])
    return parser.check_completion(failure_exit)

def track(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.annotation_track(options)

def read(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.read_or_notify()

def stop(options):
    job_pid_file = os.path.join(options['working_directory'],
                                JOB_TYPE + '.auto.pid')
    util.stop_TFLOW_process(job_pid_file, JOB_TYPE)

def clean(options):
    suffixes = ['.auto.blastx.sh', '.auto.make_db.sh']
    files = []
    out_files = [options['blast_result_file'], MATCH_PREFIX+'.All.annDB', 
                 MATCH_PREFIX+'.Best.annDB', ANNOTATION_PREFIX+'.All.annDB']
    if 'evalue_cutoffs' in options and isinstance(options['evalue_cutoffs'], list):
        for evalue_cutoff in options['evalue_cutoffs']:
            out_files.append(ANNOTATION_PREFIX + '.' + evalue_cutoff + '.annDB')
    if 'db_title' in options:
        for suffix in ['.phr', '.pin', '.psq']:
            files.append(options['db_title'] + suffix)

    remove_outfile = (options['mode'] == 'reset')
    util.clean_TFLOW_auto_files(options['job_type'], options['project_directory'],
                                options['working_directory'], remove_outfile=remove_outfile, 
                                confirm=options['confirm'], suffixes=suffixes, files=files,
                                out_files=out_files)

def test(options, silent=False):
    all_output = ''
    for job_type, command_list in [(JOB_TYPE+':BLAST', 'blast_command_list'),
                                   (JOB_TYPE+':Make_Blast_DB', 'blast_db_command_list')]:

        try:
            process = subprocess.Popen(options[command_list] + [options['test_command']], 
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            process.wait()
            output, error = process.communicate()
            all_output += output
            print ' -- %s Found!' % job_type

        except OSError as error:
            if silent:
                return False
            print ('%s Cannot Be Found ' % job_type
                   + ' With Shell Command: "%s"' % ' '.join(options[command_list]))
            if PROGRAM_URL:
                print 'If Not Installed, %s Can be Downloaded From:\n%s' % (JOB_TYPE, PROGRAM_URL)
            all_output += 'Error Number: %s\nError Text:\n%s' % (str(error.errno), error.strerror)

    return all_output

def run(options):
    if __name__ != '__main__' and options['is_pipe']:
        out_file_stream = open(options['out_file'], 'w')
        terminal_out, terminal_error = sys.stdout, sys.stderr
        sys.stdout, sys.stderr = out_file_stream, out_file_stream

    #Ensure Required Settings in Options
    for required_option in REQUIRED_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for %s not given.' % (required_option, JOB_TYPE))

    #Ensure A Type of Query Sequence File is Given
    if not any(x in options for x in ['absolute_input_analysis_file',
                                      'rel_input_analysis_file',
                                      'input_analysis_file',
                                      'result_name_file']):
        print_exit('Either input_analysis_file, absolute_input_analysis_file,'
                   + ' rel_input_analysis_file, or result_name_file paramater required.')

    #Ensure A Type of Reference Sequence File is Given
    if not any(x in options for x in ['absolute_input_reference_file',
                                      'rel_input_reference_file',
                                      'input_reference_file,'
                                      'reference_name_file']):
        print_exit('Either input_reference_file, absolute_input_reference_file,'
                   + ' rel_input_reference_file, or reference_name_file paramater required.')

    #Ensure Working Directory Exists
    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    #Assign Correct Input File Name
    if 'input_analysis_file' in options:
        if os.path.isabs(options['input_analysis_file']):
            full_input_file = options['input_analysis_file']
        else:
            full_input_file = os.path.join(options['working_directory'], 
                                           options['input_analysis_file'])
                                                   
    elif 'absolute_input_analysis_file' in options:
        full_input_file = options['absolute_input_analysis_file']

    elif 'rel_input_analysis_file' in options:
        full_input_file = os.path.join(options['project_directory'], 
                                       options['rel_input_analysis_file'])

    elif 'result_name_file' in options:
        full_result_name_file = os.path.join(options['project_directory'], 
                                             options['result_name_file'])
        if os.path.isfile(full_result_name_file):
            print ('Reading Result Sequence File Name from Provided '
                   + 'File: %s' % full_result_name_file ) 
        else:
            print_exit('Provided File: %s Containing' % full_result_name_file  
                       + ' Result Sequence File Name Not Found.')

        rf = open(full_result_name_file, 'r')
        full_input_file = rf.read().strip()
        rf.close()

        if os.path.isfile(full_input_file):
            print 'Read Result Sequence File Name: %s' % full_input_file
            print 'File Found!'
            print ''
        else:
            print_exit('Cannot Find Read Result Sequence File: %s' % full_input_file)

    input_file = os.path.basename(full_input_file)

    #Assign Correct Reference File Name
    if 'input_reference_file' in options:
        if os.path.isabs(options['input_reference_file']):
            full_reference_file = options['input_reference_file']
        else:
            full_reference_file = os.path.join(options['working_directory'], 
                                               options['input_refernece_file'])
                                                   
    elif 'absolute_input_reference_file' in options:
        full_reference_file = options['absolute_input_reference_file']

    elif 'rel_input_reference_file' in options:
        full_reference_file = os.path.join(options['project_directory'], 
                                           options['rel_input_reference_file'])

    elif 'reference_name_file' in options:
        full_reference_name_file = os.path.join(options['project_directory'], 
                                                options['reference_name_file'])
        if os.path.isfile(full_reference_name_file):
            print ('Reading Reference Sequence File Name from Provided '
                   + 'File: %s' % full_reference_name_file ) 
        else:
            print_exit('Provided File: %s Containing' % full_reference_name_file  
                       + ' Reference Sequence File Name Not Found.')

        rf = open(full_result_name_file, 'r')
        full_reference_file = rf.read().strip()
        rf.close()

        if os.path.isfile(full_reference_file):
            print 'Read Result Reference Sequence File Name: %s' % full_reference_file
            print 'File Found!'
            print ''
        else:
            print_exit('Cannot Find Read Reference Sequence File: %s' % full_reference_file)

    reference_file = os.path.basename(full_reference_file)

    #If Selected Reference File is Zipped, Unzip it
    if (full_reference_file.endswith('.gz') and os.path.isfile(full_reference_file)
        and not os.path.isfile_full_reference_file[:-3])
        print '\nSelected Reference File: %s is Zipped.' % full_reference_file 
        print 'Unzipping...'
        print ''
        sys.stdout.flush()       
        with (gzip.open(full_reference_file, 'r') as zipped_reference,
              open(full_reference_file[:-3], 'w') as unzipped_reference):
            unzipped_reference.writelines(zipped_reference)

        print ('Unzipping Complete. Setting Reference File to '
               +'Unzipped File: %s' full_reference_file[:-3])
        print ''
        full_reference_file = full_reference_file[:-3]

                                 
    #Check that Input File Exists
    if not os.path.isfile(full_input_file):
        print_exit('Input Sequence File: %s Not Found.' % full_input_file)

    #Check that Reference File Exists
    if not os.path.isfile(full_reference_file):
        print_exit('Input Reference Sequence File: %s Not Found.' % full_reference_file)

    #If Selected, Copy Input File to Working Directory
    if options['copy_input_file']:
        print ('Copying Input File: %s' % input_file
               + ' to Working Directory: %s' % options['working_directory'])
 
        working_input_file = os.path.join(options['working_directory'], input_file)
        shutil.copyfile(full_input_file, working_input_file)

        if not os.path.isfile(working_input_file):
            print_exit('Copying of File: %s to Name: %s Unsuccesful.' % (full_input_file, 
                                                                         working_input_file))
    else:
        print 'Using Input File: %s' % full_input_file 
        working_input_file = full_input_file
        print ('Input File Has %i ' % util.count_FASTA_all(working_input_file) +
               ' Detected Query Sequences.') 
               
    if options['reference_type'] in REFERENCE_TYPES:
        print 'Reference Type: %s (%s) Selected.' % (options['reference_type'], 
                                                     REFERENCE_TYPES[options['reference_type']])
        ref_type = REFERENCE_TYPES[options['reference_type']]
    else:
        print_exit(['Reference Type: %s Not Allowed.' % options['reference_type'],
                    'Options: %s' % ', '.join(REFERENCE_TYPES.keys())])

    if 'db_title' in options:
        print 'Setting Database Title to: %s' % options['db_title']
        db_title = options['db_title']
    else:
        db_title = 'BLAST_DB'

    #Prepare Blast Database
    db_command_list = options['blast_db_command_list'][:]
    db_command_list += ['-in', full_reference_file, '-dbtype', ref_type, '-title', db_title,
                        '-out', db_title]
    db_command = ' '.join(db_command_list)

    if options['write_command']:
        command_file = os.path.join(options['working_directory'],
                                    JOB_TYPE + '.auto.make_db.sh')
        util.write_file(command_file, '#!/bin/sh\n' + db_command)

    print ''
    print 'Running Command:\n    ' + db_command
    sys.stdout.flush()

    try:
        process = subprocess.Popen(db_command_list, stdout=sys.stdout, stderr=sys.stderr,
                                   cwd=options['working_directory'])
        if options['write_pid']:
            pid_file_name = os.path.join(options['working_directory'],
                                         options['job_type'] + '.auto.pid')
            util.write_file(pid_file_name, str(process.pid))
        process.wait()
        if options['write_pid']:
            util.delete_pid_file(pid_file_name)
        sys.stdout.flush()

    except KeyboardInterrupt:
        if __name__ != '__main__' and options['is_pipe']:
            sys.stdout, sys.stderr = terminal_out, terminal_error
            out_file_stream.close()

        print 'Killing %s Process.' % JOB_TYPE
        process.kill()
        raise


    print ''
    print 'Looking For Maximum Number of Annotations: %s' % options['max_matches']

    #Prepare BLAST Sequence Comparison Command
    command_list = list(options['blast_command_list'])
    command_list += ['-db', db_title, '-query', full_input_file, '-outfmt', '7', '-evalue',
                     options['blast_evalue'], '-num_threads', options['max_CPU'], 
                     '-max_target_seqs', str(options['max_matches']), '-out', 
                     options['blast_result_file']]
    command = ' '.join(command_list)

    #If Selected, Write Command to File
    if options['write_command']:
        command_file = os.path.join(options['working_directory'], JOB_TYPE + '.auto.blastx.sh')
        util.write_file(command_file, '#!/bin/sh\n' + command)

    #Perform BLAST Sequence Comparisons
    print ''
    print 'Running Command:\n    ' + command
    sys.stdout.flush()
    try:
        process = subprocess.Popen(command_list, stdout=sys.stdout, stderr=sys.stderr,
                                   cwd=options['working_directory'])
        if options['write_pid']:
            pid_file_name = os.path.join(options['working_directory'],
                                         options['job_type'] + '.auto.pid')
            util.write_file(pid_file_name, str(process.pid))
        process.wait()
        if options['write_pid']:
            util.delete_pid_file(pid_file_name)
        sys.stdout.flush()

    except KeyboardInterrupt:
        if __name__ != '__main__' and options['is_pipe']:
            sys.stdout, sys.stderr = terminal_out, terminal_error
            out_file_stream.close()

        print 'Killing %s Process.' % JOB_TYPE
        process.kill()
        raise

    print ''
    print 'Blast Completed with Out File: %s' % options['blast_result_file']
    print ''
    analyze(options)
    print ''
    print 'Annotation Complete'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout, sys.stderr = terminal_out, terminal_error
        out_file_stream.close()

#Analyze Results of Sequence Comparison
def analyze(options):
    analysis = print_return(['Performing Annotation Analysis on BLAST Result.', ''])

    #Ensure Required Settings in Options
    for required_option in REQUIRED_ANALYSIS_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for %s analysis not given.' % (required_option, 
                                                                           JOB_TYPE))
    #Ensure Working Directory Exists
    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])


    #Ensure A Type of Query Sequence File is Given
    if not any(x in options for x in ['absolute_input_analysis_file',
                                      'rel_input_analysis_file',
                                      'input_analysis_file',
                                      'result_name_file']):
        print_exit('Either input_analysis_file, absolute_input_analysis_file,'
                   + ' rel_input_analysis_file, or result_name_file paramater required.')

    #Ensure A Type of Reference Sequence File is Given
    if not any(x in options for x in ['absolute_input_reference_file',
                                      'rel_input_reference_file',
                                      'input_reference_file,'
                                      'reference_name_file']):
        print_exit('Either input_reference_file, absolute_input_reference_file,'
                   + ' rel_input_reference_file, or reference_name_file paramater required.')

    #Ensure Working Directory Exists
    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    #Assign Correct Input File Name
    if 'input_analysis_file' in options:
        if os.path.isabs(options['input_analysis_file']):
            full_input_file = options['input_analysis_file']
        else:
            full_input_file = os.path.join(options['working_directory'], 
                                           options['input_analysis_file'])
                                                   
    elif 'absolute_input_analysis_file' in options:
        full_input_file = options['absolute_input_analysis_file']

    elif 'rel_input_analysis_file' in options:
        full_input_file = os.path.join(options['project_directory'], 
                                       options['rel_input_analysis_file'])

    elif 'result_name_file' in options:
        full_result_name_file = os.path.join(options['project_directory'], 
                                             options['result_name_file'])
        if not os.path.isfile(full_result_name_file):
            print_exit('Provided File: %s Containing' % full_result_name_file  
                       + ' Result Sequence File Name Not Found.')

        rf = open(full_result_name_file, 'r')
        full_input_file = rf.read().strip()
        rf.close()

        if not os.path.isfile(full_input_file):
            print_exit('Cannot Find Read Result Sequence File: %s' % full_input_file)

    input_file = os.path.basename(full_input_file)

    #Assign Correct Reference File Name
    if 'input_reference_file' in options:
        if os.path.isabs(options['input_reference_file']):
            full_reference_file = options['input_reference_file']
        else:
            full_reference_file = os.path.join(options['working_directory'], 
                                               options['input_refernece_file'])
                                                   
    elif 'absolute_input_reference_file' in options:
        full_reference_file = options['absolute_input_reference_file']

    elif 'rel_input_reference_file' in options:
        full_reference_file = os.path.join(options['project_directory'], 
                                           options['rel_input_reference_file'])

    elif 'reference_name_file' in options:
        full_reference_name_file = os.path.join(options['project_directory'], 
                                                options['reference_name_file'])
        if not os.path.isfile(full_reference_name_file):
            print_exit('Provided File: %s Containing' % full_reference_name_file  
                       + ' Reference Sequence File Name Not Found.')

        rf = open(full_result_name_file, 'r')
        full_reference_file = rf.read().strip()
        rf.close()

        if not os.path.isfile(full_reference_file):
            print_exit('Cannot Find Read Reference Sequence File: %s' % full_reference_file)

    reference_file = os.path.basename(full_reference_file)
                                 
    #Check that Input File Exists
    if not os.path.isfile(full_input_file):
        print_exit('Input Sequence File: %s Not Found.' % full_input_file)

    #Check that Reference File Exists
    if not os.path.isfile(full_reference_file):
        print_exit('Input Reference Sequence File: %s Not Found.' % full_reference_file)

    blast_file = options['blast_result_file']
    full_blast_file = os.path.join(options['working_directory'], blast_file)

    #Check that Blast File Exists
    if not os.path.isfile(full_blast_file):
        print_exit('Input Reference Sequence File: %s Not Found.' % full_blast_file)

    if 'name_map_file' in options and options['name_map_file']:
        if os.path.isabs(options['name_map_file']):
            full_name_map_file = options['name_map_file']
            name_map_file = os.path.basename(full_name_map_file)
        else:
            name_map_file = options['name_map_file']
            full_name_map_file = os.path.join(options['working_directory'], name_map_file)
        
        analysis += print_return(['Reading Name Map File:', '    %s' % full_name_map_file])
    else: 
        full_name_map_file = None
        name_map_file = None

    #If Provided, Check that Name Map File Exists
    if name_map_file:
        if not os.path.isfile(full_blast_file):
            print_exit('Provided Name Mapping File: %s Not Found.' % full_name_map_file)

    db = Annotation_Database()
    analysis += print_return(['Beginning Annotation...', ''])

    #Read # of Sequences in Input File
    input_sequence_count = util.count_FASTA_all(full_input_file)

    analysis += print_return(['Total Sequences in input file %s:' % full_input_file
                              + ' %i ' % input_sequence_count, ''])

    #Read # of Lines in File
    with open(full_blast_file) as blast_file_object:
        for total_lines, line in enumerate(blast_file_object, start=1):
            pass

    analysis += print_return(['Total Lines in file %s: %s ' % (full_blast_file, str(total_lines)),
                              ''])

    #Read Blast File Outputs and Count Genes Found Over Threshold
    blast_file_object = open(full_blast_file, 'r')

    print_counter = 0
    NUM_PRINTS = 1000
    print_counter_threshold = total_lines/NUM_PRINTS

    db_len = 0
    for (line_number, line) in enumerate(blast_file_object, start=1):
        print_counter += 1
        if line.startswith('#'):
            continue
        split_line = line.split()
        if not split_line:
            print_exit('Blank Line Found in Blast Results File at Line Number %i' % line_number)
        elif len(split_line) < 11:
            print_exit([('Problem with formatting of line number %i ' % line_number
                        + 'in blast results file: %s' % full_blast), 'Line:', line.strip()])
        query_sequence = split_line[0]
        match_sequence = split_line[1]
        e_score_string = split_line[10]
        e_score = float(e_score_string)
        record = Annotation_Record(ID=match_sequence, eVal=e_score, fileName=reference_file, 
                                   annotation=None)
        db.add_record(query_sequence, record)

        if print_counter >= print_counter_threshold:
            db_len = len(db)
            print ('\r%s Matched Sequences Found,' % (str(db_len))
                   + ' %s complete.          ' % util.percent_string(line_number, total_lines)), 
            sys.stdout.flush()
            print_counter = (print_counter_threshold - print_counter)

    print '\r' + ' ' * 79 + '\r',
    analysis += print_return('%s Matched Sequences Found.' % (str(db_len)))

    blast_file_object.close()
    if options['write_all_matches']:       
        analysis += print_return('Writing All Sequence Matches to File: ' + 
                                 MATCH_PREFIX + '.All.annDB')
        sys.stdout.flush()
        db.write(MATCH_PREFIX + '.All.annDB')
        annotation_count, record_count = db.count()
        analysis += print_return(['%i Sequences ' % annotation_count
                                  + ' with %i Total Matches Written.' % record_count, ''])

    db.cull(subset='best')
    if options['write_best_matches']:       
        analysis += print_return('Writing Best Sequence Matches to File: ' + 
                                 MATCH_PREFIX + '.Best.annDB')
        sys.stdout.flush()
        db.write(MATCH_PREFIX + '.Best.annDB')
        annotation_count, record_count = db.count()
        analysis += print_return(['%i Sequences ' % annotation_count
                                  + ' with %i Total Matches Written.' % record_count, ''])

    if name_map_file:
        input_sequence_count = 0
        analysis += print_return('Reading Provided Name Map: %s' % full_name_map)
        name_map = Name_Map(full_name_map)
        analysis += print_return('Remapping Sequence Names to Name Map...')
        db.map_names(name_map, debug=False)      
        annotation_count, record_count = db.count()
        analysis += print_return(['%i Sequences Remapped ' % annotation_count 
                                 + 'with %i Matches.'% record_count, ''])

        if options['write_best_matches']:       
            analysis += print_return('Writing Best Matches to File: ' + 
                                     MATCH_PREFIX + '.Remapped.Best.annDB')
            sys.stdout.flush()
            db.write(MATCH_PREFIX + '.Remapped.Best.annDB')
            annotation_count, record_count = db.count()
            analysis += print_return(['%i Sequences ' % annotation_count
                                      + ' with %i Total Matches Written.' % record_count, ''])

    analysis += print_return('Reading Annotation Strings from Reference File: %s' % full_reference_file)
    annotation_map = Annotation_Map(full_reference_file)

    analysis += print_return('Mapping Annotations to Matches...')
    map_counts = db.map_annotations(annotation_map)
    analysis += print_return(['%i Total Annotations Mapped.' % map_counts, ''])

    analysis += print_return('Writing All Annotations to File: ' + 
                             ANNOTATION_PREFIX + '.All.annDB')
    sys.stdout.flush()
    db.write(ANNOTATION_PREFIX + '.All.annDB')
    annotation_count, record_count = db.count()
    analysis += print_return(['%i Sequences ' % annotation_count
                              + ' with %i Total Matches Written.' % record_count, ''])


    thresholds = options['evalue_cutoffs']
    last_threshold = 10
    for threshold in thresholds:
        if float(threshold) >= float(last_threshold):
            print_exit(['Thresholds: %s ' %', '.join(thresholds)
                        +'Must Be in Descending Order.',
                        '(1e-10, then 1e-20, then 1e-40, etc.)'])
        else:
            last_threshold = threshold

    analysis += print_return(['Analyzing Annotations for Evalue Thresholds: '
                              + ', '.join(options['evalue_cutoffs']), ''])

    report_dicts = []
    formatted_reports = []
    for threshold in thresholds:
        (ini_seqs, ini_records, final_seqs, final_records) = db.cull(threshold=threshold)

        analysis += print_return('Writing Theshold %s Annotations to File: ' % threshold + 
                                 ANNOTATION_PREFIX + '.' + threshold + '.annDB')
        sys.stdout.flush()
        db.write(ANNOTATION_PREFIX + '.' + threshold + '.annDB')
        annotation_count, record_count = db.count()
        analysis += print_return(['%i Sequences' % annotation_count
                                  + ' with %i Total Matches Written.' % record_count, ''])

        if input_sequence_count:
            formatted_input_sequence_count = str(input_sequence_count)
            percent = util.percent_string(final_seqs, input_sequence_count)
        else:
            formatted_input_sequence_count = '-'
            percent = 'N/A%'

        headers = ['Analys.', 'Cutoff', 'TotSqs.', 'AnnSqs.', 'Percent', 'TotAnn.', 'Remapd.']
        data_grid = ['Annot.', threshold, formatted_input_sequence_count, final_seqs, 
                     percent, final_records, bool(name_map_file)]
        formatted_data = [str(x) for x in data_grid]
    
        report_dict = dict(zip(headers, formatted_data))
        report_dict['report_type'] = 'annotation'
        report_dicts.append(report_dict)
        formatted_reports.append(formatted_data)

    analysis += print_return(['', 'Tab Separated Output:', '\t'.join(headers)])
    for formatted_report in formatted_reports:
        analysis += print_return('\t'.join(formatted_report))
    analysis += print_return('')
       
    #If Selected, Write Analysis Report
    if options['write_report']:
        report_file = os.path.join(options['working_directory'],
                                   JOB_TYPE + '.report')
        if len(report_dicts) > 1:
            aux_reports = report_dicts[1:]
        else:
            aux_reports = []
        util.write_report(report_file, report_dict, aux_reports=aux_reports)

    return analysis

            
    
    
    

