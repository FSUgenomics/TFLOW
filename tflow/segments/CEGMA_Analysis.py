#TFLOW Segment: Analyze FASTA File for Gene Recapture using CEGMA Benchmark Database
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

if __name__ == "__main__" or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../'))
    import tflow.segments
    __package__ = "tflow.segments"

from .. import local_settings
from .parser_class import OutputParser
from ..util import (print_exit, print_except, write_file, write_report, delete_pid_file, 
                    percent_string, stop_TFLOW_process)
from .. import util

if hasattr(local_settings, 'CEGMA_FILE'):
    CEGMA_FILE = local_settings.CEGMA_FILE
else:
    CEGMA_FILE = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..',
                              'sequence_files', 'CEGMA.fasta') 

if hasattr(local_settings, 'BLAST_LOCATION'):
    BLAST_LOCATION = local_settings.BLAST_LOCATION
else:
    BLAST_LOCATION = ''

if hasattr(local_settings, 'MAKE_BLAST_DB_LOCATION'):
    MAKE_BLAST_DB_LOCATION = local_settings.MAKE_BLAST_DB_LOCATION
else:
    MAKE_BLAST_DB_LOCATION = ''

JOB_TYPE = 'CEGMA_Analysis'
PROGRAM_URL = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download'
DATASET_DETAILS ='''Core Eukaryotic Genes Mapping Approach (CEGMA): Core Eukaryotic Genes Dataset
Version:  Acquired 2015-04-22
URL:      http://korflab.ucdavis.edu/datasets/cegma/
Citation: Genis Parra, Keith Bradnam and Ian Korf. CEGMA: a pipeline to accurately annotate
          core genes in eukaryotic genomes." Bioinformatics, 23: 1061-1067 (2007)
'''
SEGMENT_FOR_VERSION = '2.2.29'
BLAST_COMMAND = os.path.join(BLAST_LOCATION, 'blastx')
BLAST_COMMAND_LIST = [BLAST_COMMAND]
BLAST_DB_COMMAND = os.path.join(MAKE_BLAST_DB_LOCATION, 'makeblastdb')
BLAST_DB_COMMAND_LIST = [BLAST_DB_COMMAND]
TEST_COMMAND = '-h'
OUT_FILE = 'CEGMA_Analysis.out'
MILESTONES = ['CEGMA Benchmarking Analysis Complete']
TERMINAL_FLAGS = ['CEGMA Analysis Done']
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Exception: ERROR',
                 'Not Found']
DEFAULT_SETTINGS = {'working_directory':'CEGMA_Analysis',
                    'CEGMA_file':CEGMA_FILE,
                    'copy_input_file':True,
                    'max_CPU':'4',
                    'evalue':'1e-5',
                    'evalue_cutoff': '1e-20',
                    'blast_result_file':'blast.out',
                    'print_missing_genes':False,
                    'print_matches':False,
                    #TFLOW CEGMA_Analysis Settings
                    'blast_command':BLAST_COMMAND,
                    'blast_command_list':BLAST_COMMAND_LIST,
                    'blast_db_command':BLAST_DB_COMMAND,
                    'blast_db_command_list':BLAST_DB_COMMAND_LIST,
                    'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    'dataset_details':DATASET_DETAILS,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    'write_report':True,
                    'write_command':True,
                    'write_pid':True,
                    }

REQUIRED_SETTINGS = ['blast_command_list', 'blast_db_command_list', 'working_directory', 
                     'CEGMA_file', 'copy_input_file', 'evalue', 'max_CPU', 'blast_result_file',
                     'evalue_cutoff', 'print_missing_genes', 'write_command', 'write_report', 
                     'write_pid', 'print_matches']

REQUIRED_ANALYSIS_SETTINGS = ['CEGMA_file', 'blast_result_file', 'evalue_cutoff', 
                              'working_directory', 'print_missing_genes', 'write_report', 
                              'print_matches']

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
    job_pid_file = os.path.join(options['working_directory'],
                                JOB_TYPE + '.auto.pid')
    stop_TFLOW_process(job_pid_file, JOB_TYPE)

def clean(options):
    files = ['CEGMA.phr', 'CEGMA.psq', 'CEGMA_Make_DB.auto.sh', 
             'CEGMA.pin', 'CEGMA_tblastn.auto.sh']
    out_files = [options['blast_result_file']]
    remove_outfile = (options['mode'] == 'reset')
    util.clean_TFLOW_auto_files(options['job_type'], options['project_directory'],
                                options['working_directory'], remove_outfile=remove_outfile, 
                                confirm=options['confirm'], files=files, out_files=out_files)

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

    #Ensure A Type of Input File is Given
    if not any(x in options for x in ['absolute_input_analysis_file',
                                      'rel_input_analysis_file',
                                      'result_name_file']):
        print_exit('Either absolute_input_analysis_file, rel_input_analysis_file, or'
                   + ' result_name_file paramater required.')

    #Ensure Working Directory Exists
    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    #Assign Correct Input File Name
    if 'absolute_input_analysis_file' in options:
        full_input_file = options['absolute_input_analysis_file']
        input_file = os.path.basename(full_input_file)

    elif 'rel_input_analysis_file' in options:
        full_input_file = os.path.join(options['project_directory'], 
                                       options['rel_input_analysis_file'])
        input_file = os.path.basename(options['rel_input_analysis_file'])

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

    #Ensure Input File Exists
    if not os.path.isfile(full_input_file):
        print_exit('Input File: %s Not Found.' % full_input_file)

    #Print Dataset Details
    if 'dataset_details' in options:
        print 'Details on Benchmarking Dataset:'
        print options['dataset_details']
        print ''

    #If CEGMA Sequence File is Zipped, Unzip it
    if os.path.isfile(options['CEGMA_file'] +'.gz'):
        print ('\nProvided CEGMA File: %s' % options['CEGMA_file']
               + 'Found in Zipped Format: %s' % options['CEGMA_file'] + '.gz')
        print 'Unzipping...'
        print ''
        process = subprocess.Popen(['gunzip', options['CEGMA_file'] +'.gz'], stdout=sys.stdout,
                                   stderr=sys.stderr, cwd=options['working_directory'])
        if options['write_pid']:
            pid_file_name = os.path.join(options['working_directory'],
                                         options['job_type'] + '.auto.pid')
            write_file(pid_file_name, str(process.pid))

        process.wait()

        if options['write_pid']:
            delete_pid_file(pid_file_name)

        sys.stdout.flush()
        print ''

    #Ensure CEGMA Sequence File Exists
    if os.path.isfile(options['CEGMA_file']):
        print 'Provided CEGMA Sequence File Found: %s' % options['CEGMA_file']
    else:
        print_exit('Provided CEGMA Sequence File: %s Cannot Be Found.' % options['CEGMA_file'])

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


    #Prepare Blast Database
    db_command_list = list(options['blast_db_command_list'])
    db_command_list += ['-in', options['CEGMA_file'], '-dbtype', 'prot', '-title', 'CEGMA',
                     '-out', 'CEGMA']
    db_command = ' '.join(db_command_list)

    if options['write_command']:
        command_file = os.path.join(options['working_directory'],
                                    'CEGMA_Make_DB.auto.sh')
        write_file(command_file, '#!/bin/sh\n' + db_command)

    print ''
    print 'Running Command:\n    ' + db_command

    sys.stdout.flush()

    try:
        process = subprocess.Popen(db_command_list, stdout=sys.stdout, stderr=sys.stderr,
                                   cwd=options['working_directory'])
        if options['write_pid']:
            pid_file_name = os.path.join(options['working_directory'],
                                         options['job_type'] + '.auto.pid')
            write_file(pid_file_name, str(process.pid))

        process.wait()

        if options['write_pid']:
            delete_pid_file(pid_file_name)

        sys.stdout.flush()

    except KeyboardInterrupt:
        if __name__ != '__main__' and options['is_pipe']:
            sys.stdout, sys.stderr = terminal_out, terminal_error
            out_file_stream.close()

        print 'Killing %s Process.' % JOB_TYPE
        process.kill()
        raise

    #Prepare BLAST Sequence Comparison Command
    command_list = list(options['blast_command_list'])
    command_list += ['-db', 'CEGMA', '-query', full_input_file, '-outfmt', '6', '-evalue',
                     options['evalue'], '-num_threads', options['max_CPU'], '-out', 
                     options['blast_result_file']]
    command = ' '.join(command_list)

    #If Selected, Write Command to File
    if options['write_command']:
        command_file = os.path.join(options['working_directory'], 'CEGMA_blastx.auto.sh')
        write_file(command_file, '#!/bin/sh\n' + command)

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
            write_file(pid_file_name, str(process.pid))

        process.wait()

        if options['write_pid']:
            delete_pid_file(pid_file_name)

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
    print 'CEGMA Benchmarking Analysis Complete'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout, sys.stderr = terminal_out, terminal_error
        out_file_stream.close()

#Analyze Results of Sequence Comparison
def analyze(options):
    analysis = 'Analyzing CEMGA Recapture BLAST Result.\n\n'

    #Ensure Required Settings in Options
    for required_option in REQUIRED_ANALYSIS_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for %s not given.' % (required_option, JOB_TYPE))

    #Ensure Working Directory Exists
    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    #Prepare File Names
    full_blast = os.path.join(options['working_directory'], options['blast_result_file'])
    full_cegma = options['CEGMA_file']

    #Ensure Necessary Files Exists
    if not os.path.isfile(full_blast):
        print_exit('Blast Result File: %s Not Found.' % full_blast)

    if not os.path.isfile(full_cegma):
        print_exit('CEGMA Sequence File: %s Not Found.' % full_cegma)

    analysis = '\nAnalyzing Blast Result File %s\n' % full_blast
    analysis += '    With CEGMA file: %s\n' % full_cegma

    #Read Expected Genes
    genes = {}
    cegma_file = open(full_cegma, 'r')
    for line in cegma_file:
        if line.startswith('>'):
            gene = line.split('___')[-1].strip()
            genes[gene]=False
    cegma_file.close()

    expected_gene_count = len(genes)
    analysis += '\nExpected Genes: %i\n' % expected_gene_count

    cutoff_float = float(options['evalue_cutoff'])

    #Read Blast File Outputs and Count Genes Found Over Threshold
    blast_file = open(full_blast, 'r')
    for (line_number, line) in enumerate(blast_file,start=1):
        split_line = line.split()
        if not split_line:
            print_exit('Blank Line Found in Blast Results File at Line Number %i' % line_number)
        elif len(split_line) < 11:
            print_exit([('Problem with formatting of line number %i ' % line_number
                        + 'in blast results file: %s' % full_blast), 'Line:', line.strip()])
        sequence = split_line[0]
        gene = split_line[1].split('___')[-1].strip()
        e_score_string = split_line[10]
        e_score = float(e_score_string)

        #Mark Gene as Present if Hit Exists over Threshold Value
        if e_score <= cutoff_float:
            if options['print_matches'] and not genes[gene]:
                analysis += 'Match: %s %s %s\n' % (sequence, gene, e_score_string)
            genes[gene] = True

    #Count Number of Found and Missing Genes
    found_gene_count = 0
    missing_genes = []
    for gene in genes:
        if genes[gene]:
            found_gene_count += 1
        else:
            missing_genes.append(gene)
        missing_gene_count = len(missing_genes)

    #Ensure that Found/Missing Genes Sums to Expected Total
    if missing_gene_count + found_gene_count != expected_gene_count:
        print_except('PROBLEM!, Found: %i + ' % found_gene_count
                     + 'Missing: %i Genes != Expected: %i' % (missing_gene_count,
                                                              expected_gene_count)) 
                                                                                  
    #Report Results
    analysis += 'Genes Found: %i\n' % found_gene_count
    analysis += 'Genes Missing: %i\n' % missing_gene_count
    if options['print_missing_genes'] and missing_genes:
        analysis += 'Missing Genes: ' + ' '.join(missing_genes) + '\n'

    percent = percent_string(found_gene_count, expected_gene_count)
    analysis += 'Percent CEGMA Genes Present: %s\n' % percent

    headers = ['Analys.', 'Cutoff', 'Expect.', 'Found', 'Missing', 'Total', 'Percent']
    data_grid = ['CEGMA', options['evalue_cutoff'], expected_gene_count, found_gene_count, 
                 missing_gene_count, expected_gene_count, percent]
    formatted_data = [str(x) for x in data_grid]

    analysis += '\n'
    analysis += 'Tab Separated Output:\n'
    analysis += '\t'.join(headers) + '\n'
    analysis += '\t'.join(formatted_data) + '\n'

    report_dict = dict(zip(headers, formatted_data))
    report_dict['report_type'] = 'recapture'

    #If Selected, Write Analysis Report
    if options['write_report']:
        report_file = os.path.join(options['working_directory'],
                                   JOB_TYPE + '.report')
        write_report(report_file, report_dict)

    print analysis
    return analysis

            
    
    
    

