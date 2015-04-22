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
    import tflow.segments
    __package__ = "tflow.segments"

from .. import local_settings
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

from .parser_class import OutputParser
from ..util import print_exit, print_except, write_file, percent_string

JOB_TYPE = 'CEGMA_Analysis'
PROGRAM_URL = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download'
SEGMENT_FOR_VERSION = '2.2.29'
BLAST_COMMAND = os.path.join(BLAST_LOCATION, 'blastx')
BLAST_COMMAND_LIST = [BLAST_COMMAND]
BLAST_DB_COMMAND = os.path.join(MAKE_BLAST_DB_LOCATION, 'makeblastdb')
BLAST_DB_COMMAND_LIST = [BLAST_DB_COMMAND]
TEST_COMMAND = '-h'
OUT_FILE = 'CEGMA_Analysis.out'
MILESTONES = ['CEGMA Analysis Done']
TERMINAL_FLAGS = []
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Not Found']
DEFAULT_SETTINGS = {'working_directory':'CEGMA_Analysis',
                    'CEGMA_file':CEGMA_FILE,
                    'copy_input_file':True,
                    'max_CPU':'4',
                    'evalue':'1e-5',
                    'evalue_cutoff':'1e-20',
                    'blast_result_file':'blast.out',
                    'print_missing_genes':False,
                    #TFLOW CEGMA_Analysis Settings
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
                    }

REQUIRED_SETTINGS = ['blast_command_list', 'blast_db_command_list', 'working_directory', 
                     'CEGMA_file', 'copy_input_file', 'evalue', 'max_CPU', 'blast_result_file',
                     'evalue_cutoff', 'print_missing_genes', 'write_command', 'write_report']

REQUIRED_ANALYSIS_SETTINGS = ['CEGMA_file', 'blast_result_file', 'evalue_cutoff', 
                              'print_missing_genes', 'write_report' ]

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

def read(options):
    parser = Parser()
    parser.out_file = options['out_file']
    parser.read_or_notify()

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

    for required_option in REQUIRED_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for %s not given.' % (required_option, JOB_TYPE))

    if not any(x in options for x in ['absolute_input_analysis_file',
                                      'rel_input_analysis_file']):
        print_exit('Either absolute_input_analysis_file or rel_input_analysis_file'
                   + ' paramater required.')

    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    if 'absolute_input_analysis_file' in options:
        full_input_file = options['absolute_input_analysis_file']
        input_file = os.path.basename(full_input_file)

    elif 'rel_input_analysis_file' in options:
        full_input_file = os.path.join(options['project_directory'], 
                                       options['rel_input_analysis_file'])
        input_file = os.path.basename(options['rel_input_analysis_file'])

    if not os.path.isfile(full_input_file):
        print_exit('Input File: %s Not Found.' % full_input_file)

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

    db_command_list = options['blast_db_command_list'][:]
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
        process.wait()
        sys.stdout.flush()

    except KeyboardInterrupt:
        if __name__ != '__main__' and options['is_pipe']:
            sys.stdout, sys.stderr = terminal_out, terminal_error
            out_file_stream.close()

        print 'Killing %s Process.' % JOB_TYPE
        process.kill()
        raise


    #Perform BLAST

    command_list = options['blast_command_list'][:]
    command_list += ['-db', 'CEGMA', '-query', full_input_file, '-outfmt', '6', '-evalue',
                     options['evalue'], '-num_threads', options['max_CPU'], '-out', 
                     options['blast_result_file']]
    command = ' '.join(command_list)

    if options['write_command']:
        command_file = os.path.join(options['working_directory'], 'CEGMA_tblastn.auto.sh')
        write_file(command_file, '#!/bin/sh\n' + command)

    print ''
    print 'Running Command:\n    ' + command

    sys.stdout.flush()

    try:
        process = subprocess.Popen(command_list, stdout=sys.stdout, stderr=sys.stderr,
                                   cwd=options['working_directory'])
        process.wait()
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
    print 'CEGMA Analysis Done'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout = terminal_out
        sys.stderr = terminal_error
        out_file_stream.close()

def analyze(options):
    analysis = 'Analyzing CEMGA Recapture BLAST Result.\n\n'

    for required_option in REQUIRED_ANALYSIS_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for %s not given.' % (required_option, JOB_TYPE))

    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    full_blast = os.path.join(options['working_directory'], options['blast_result_file'])
    full_cegma = options['CEGMA_file']

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

    blast_file = open(full_blast, 'r')
    for line in blast_file:
        split_line = line.split()
        sequence = split_line[0]
        gene = split_line[1].split('___')[-1].strip()
        e_score = float(split_line[10])

        if e_score <= cutoff_float:
            genes[gene] = True

    found_gene_count = 0
    missing_genes = []
    for gene in genes:
        if genes[gene]:
            found_gene_count += 1
        else:
            missing_genes.append(gene)
    
    missing_gene_count = len(missing_genes)

    if missing_gene_count + found_gene_count != expected_gene_count:
        print_except('PROBLEM!, Found: %i + ' % found_gene_count
                     + 'Missing: %i Genes != Expected: %i' % (missing_gene_count,
                                                              expected_gene_count)) 
                                                                                  
    analysis += 'Genes Found: %i\n' % found_gene_count
    analysis += 'Genes Missing: %i\n' % missing_gene_count
    if options['print_missing_genes'] and missing_genes:
        analysis += 'Missing Genes: ' + ' '.join(missing_genes) + '\n'

    percent = percent_string(found_gene_count, expected_gene_count)
    analysis += 'Percent CEGMA Genes Present: %s\n' % percent



    headers = ['Analys.', 'Cutoff', 'Expect.', 'Found', 'Missing', 'Total', 'Percent']
    #formatted_header = [str(x).ljust(int(options['analysis_column_size'])) for x in headers]
    formatted_header = [str(x) for x in headers]

    analysis += '\n'
    analysis += 'Tab Separated Output:\n'

    report = '\t'.join(formatted_header) + '\n'

    data_grid = ['CEGMA', options['evalue_cutoff'], expected_gene_count, found_gene_count, 
                 missing_gene_count, expected_gene_count, percent]
    #formatted_data = [str(x).ljust(int(options['analysis_column_size'])) for x in data_grid]
    formatted_data = [str(x) for x in data_grid]
    report += '\t'.join(formatted_data) + '\n'

    analysis += report

    if options['write_report']:
        report_file = os.path.join(options['working_directory'],
                                   JOB_TYPE + '.report')
        write_file(report_file, report)

    print analysis
    return analysis

            
    
    
    

