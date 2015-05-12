#TFLOW Segment: Analyze FASTA File for Gene Recapture using BUSCO Benchmark Database
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

BUSCO_FILES = {'arthropoda':'BUSCO_Arthropoda.fas', 
               'vertebrata':'BUSCO_Vertebrata.fas',
               'fungi':'BUSCO_Fungi.fas',
               'metazoa':'BUSCO_Metazoa.fas'}

if __name__ == "__main__" or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../'))
    import tflow.segments
    __package__ = "tflow.segments"

from .parser_class import OutputParser
from ..util import print_exit, print_except, write_file, write_report, percent_string, lowercase
from .. import local_settings

if hasattr(local_settings, 'BUSCO_LOCATION'):
    BUSCO_LOCATION = local_settings.BUSCO_LOCATION
else:
    BUSCO_LOCATION = os.path.join(os.path.dirname(os.path.realpath(__file__)), '..', 
                                  'sequence_files')

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



JOB_TYPE = 'BUSCO_Analysis'
PROGRAM_URL = 'http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download'
DATASET_DETAILS = ''' Benchmarking sets of Universal Single-Copy Orthologs (BUSCO): Metazoa, Arthropoda,
    Vertebrata, and Fungi Datasets
Version:  OrthoDB7, Acquired 2015-04-22
URL:      ftp://cegg.unige.ch/OrthoDB7/BUSCO/
Citation: Waterhouse et al, Nucleic Acids Research, January 2013, PMID:23180791
          OrthoDB: a hierarchical catalog of animal, fungal and bacterial orthologs.
'''
SEGMENT_FOR_VERSION = '2.2.29'
BLAST_COMMAND = BLAST_EXEC
BLAST_COMMAND_LIST = [BLAST_COMMAND]
BLAST_DB_COMMAND = MAKE_BLAST_DB_EXEC
BLAST_DB_COMMAND_LIST = [BLAST_DB_COMMAND]
TEST_COMMAND = '-h'
OUT_FILE = 'BUSCO_Analysis.out'
MILESTONES = ['BUSCO Benchmarking Analysis Complete']
TERMINAL_FLAGS = ['BUSCO Analysis Done']
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Not Found',
                 ]
DEFAULT_SETTINGS = {'working_directory':'BUSCO_Analysis',
                    'BUSCO_type':'vertebrata',
                    'BUSCO_location':BUSCO_LOCATION,
                    'copy_input_file':True,
                    'max_CPU':'4',
                    'evalue':'1e-5',
                    'evalue_cutoff':'1e-20',
                    'blast_result_file':'blast.out',
                    'print_missing_genes':False,
                    'print_matches':False,
                    #TFLOW BUSCO_Analysis Settings
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
                    }

REQUIRED_SETTINGS = ['blast_command_list', 'blast_db_command_list', 'working_directory', 
                     'copy_input_file', 'evalue', 'max_CPU', 'blast_result_file', 'evalue_cutoff',
                     'print_missing_genes', 'write_command', 'write_report', 'print_matches']

REQUIRED_ANALYSIS_SETTINGS = ['blast_result_file', 'evalue_cutoff', 'print_missing_genes',
                              'print_matches', 'write_report']

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

    #Ensure a BUSCO file or type is given
    if not any(x in options for x in ['BUSCO_file', 'BUSCO_type']):
        print_exit('Either BUSCO_file or BUSCO_type paramater required.')

    #Ensure Working Directory Exists
    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    #Print Dataset Details
    if 'dataset_details' in options:
        print 'Details on Benchmarking Dataset:'
        print options['dataset_details']
        print ''

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

    #Find/Validate BUSCO File Selection
    if 'BUSCO_file' in options:
        full_BUSCO_file_name = options['BUSCO_file']
        print 'BUSCO File: %s Given.' % options['BUSCO_file']

    elif 'BUSCO_type' in options:
        if 'BUSCO_location' not in options or not options['BUSCO_location']:
            print_exit('BUSCO_type: %s Given ' % options['BUSCO_type']
                       + 'but BUSCO_location not given.')
        
        if not os.path.isdir(options['BUSCO_location']):
            print_exit('BUSCO File Location: %s Not Found.' % options['BUSCO_location'])

        BUSCO_type = lowercase(options['BUSCO_type'])
        if BUSCO_type in BUSCO_FILES:
            print 'BUSCO File Type: %s Provided.' % BUSCO_type
            full_BUSCO_file_name = os.path.join(options['BUSCO_location'], 
                                                BUSCO_FILES[BUSCO_type])
        else:
            print_exit([('Selected BUSCO Type: %s Not Available.' % BUSCO_Type),
                       'Please Select from Types:', ', '.join(BUSCO_FILES.keys())]) 

    #If Selected BUSCO File is Zipped, Unzip it
    if os.path.isfile(full_BUSCO_file_name +'.gz'):
        print ('\nSelected BUSCO File: %s' % full_BUSCO_file_name
               + 'Found in Zipped Format: %s' % full_BUSCO_file_name + '.gz')
        print 'Unzipping...'
        print ''
        sys.stdout.flush()
        process = subprocess.Popen(['gunzip', full_BUSCO_file_name +'.gz'], stdout=sys.stdout, 
                                   stderr=sys.stderr, cwd=options['working_directory'])
        process.wait()
        sys.stdout.flush()
        print ''

    #Ensure Provided/Selected BUSCO File Exists
    if os.path.isfile(full_BUSCO_file_name):
        print 'Selected BUSCO File Found: %s' % full_BUSCO_file_name
        if 'BUSCO_file' not in options:
            options['BUSCO_file'] = full_BUSCO_file_name

    else:
        print_exit('Selected BUSCO File: %s Cannot Be Found.' % full_BUSCO_file_name)
                                 
    #Check that Input File Exists
    if not os.path.isfile(full_input_file):
        print_exit('Input File: %s Not Found.' % full_input_file)

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


    #Prepare Blast Database Name
    if 'BUSCO_type' in options:
        title='BUSCO_' + options['BUSCO_type'].title()
    else:
        BUSCO_file_name = os.path.basename(options['BUSCO_file'])
        if BUSCO_file_name in BUSCO_FILES.values():
            for name in BUSCO_FILES:
                if BUSCO_file_name == BUSCO_FILES[name]:
                    title = 'BUSCO_' + name.title()
                    break
        else:
            title = 'BUSCO'
                    
    #Prepare Blast Database
    db_command_list = options['blast_db_command_list'][:]
    db_command_list += ['-in', full_BUSCO_file_name, '-dbtype', 'prot', '-title', title,
                        '-out', title]
    db_command = ' '.join(db_command_list)

    if options['write_command']:
        command_file = os.path.join(options['working_directory'],
                                    'BUSCO_Make_DB.auto.sh')
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

    #Prepare BLAST Sequence Comparison Command
    command_list = list(options['blast_command_list'])
    command_list += ['-db', title, '-query', full_input_file, '-outfmt', '6', '-evalue',
                     options['evalue'], '-num_threads', options['max_CPU'], '-out', 
                     options['blast_result_file']]
    command = ' '.join(command_list)

    #If Selected, Write Command to File
    if options['write_command']:
        command_file = os.path.join(options['working_directory'], 'BUSCO_tblastn.auto.sh')
        write_file(command_file, '#!/bin/sh\n' + command)

    #Perform BLAST Sequence Comparisons
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
    print 'BUSCO Benchmarking Analysis Complete'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout, sys.stderr = terminal_out, terminal_error
        out_file_stream.close()

#Analyze Results of Sequence Comparison
def analyze(options):
    analysis = 'Analyzing BUSCO Recapture BLAST Result.\n\n'

    #Ensure Required Settings in Options
    for required_option in REQUIRED_ANALYSIS_SETTINGS:
        if required_option not in options:
            print_exit('Required Option: %s for %s not given.' % (required_option, JOB_TYPE))

    #Ensure a BUSCO file or type is given
    if not any(x in options for x in ['BUSCO_file', 'BUSCO_type']):
        print_exit('Either BUSCO_file or BUSCO_type paramater required.')

    #Ensure Working Directory Exists
    if not os.path.isdir(options['working_directory']):
        print_exit('Working Directory: %s Not Found.' % options['working_directory'])

    #Find/Validate BUSCO File Selection
    if 'BUSCO_file' in options:
        full_BUSCO_file_name = options['BUSCO_file']
        print 'BUSCO File: %s Given.' % options['BUSCO_file']

    elif 'BUSCO_type' in options:
        if 'BUSCO_location' not in options or not options['BUSCO_location']:
            print_exit('BUSCO_type: %s Given ' % options['BUSCO_type']
                       + 'but BUSCO_location not given.')
        
        if not os.path.isdir(options['BUSCO_location']):
            print_exit('BUSCO File Location: %s Not Found.' % options['BUSCO_location'])

        BUSCO_type = lowercase(options['BUSCO_type'])
        if BUSCO_type in BUSCO_FILES:
            print 'BUSCO File Type: %s Provided.' % BUSCO_type
            full_BUSCO_file_name = os.path.join(options['BUSCO_location'], 
                                                BUSCO_FILES[BUSCO_type])
        else:
            print_exit([('Selected BUSCO Type: %s Not Available.' % BUSCO_Type),
                       'Please Select from Types:', ', '.join(BUSCO_FILES.keys())]) 

    #Ensure Provided/Selected BUSCO File Exists
    if os.path.isfile(full_BUSCO_file_name):
        print 'Selected BUSCO File Found: %s' % full_BUSCO_file_name
    else:
        print_exit('Selected BUSCO File: %s Cannot Be Found.' % full_BUSCO_file_name)
                                 
    full_blast = os.path.join(options['working_directory'], options['blast_result_file'])

    #Ensure Blast Output File Exists
    if not os.path.isfile(full_blast):
        print_exit('Blast Result File: %s Not Found.' % full_blast)

    analysis = '\nAnalyzing Blast Result File %s\n' % full_blast
    analysis += '    With BUSCO file: %s\n' % full_BUSCO_file_name

    #Read Expected Genes
    BUSCO_sequences = {}
    genes = {}
    BUSCO_file = open(full_BUSCO_file_name, 'r')
    for line in BUSCO_file:
        if line.startswith('>'):
            split_line = line.lstrip('>').split()
            sequence = split_line[0]
            gene = split_line[1]
            BUSCO_sequences[sequence] = gene
            genes[gene] = False
    BUSCO_file.close()

    expected_gene_count = len(genes)
    analysis += '\nExpected Genes: %i\n' % expected_gene_count

    cutoff_float = float(options['evalue_cutoff'])

    #Read Blast File Outputs and Count Genes Found Over Threshold
    blast_file = open(full_blast, 'r')
    for (line_number, line) in enumerate(blast_file, start=1):
        split_line = line.split()
        if not split_line:
            print_exit('Blank Line Found in Blast Results File at Line Number %i' % line_number)
        elif len(split_line) < 11:
            print_exit([('Problem with formatting of line number %i ' % line_number
                        + 'in blast results file: %s' % full_blast), 'Line:', line.strip()])
        sequence = split_line[0]
        BUSCO_sequence = split_line[1]
        if BUSCO_sequence in BUSCO_sequences:
            gene = BUSCO_sequences[BUSCO_sequence]
        else:
            print_except(['Unexpected BUSCO Sequence Hit: %s Found.' % BUSCO_sequence,
                          'Cannot Identify Gene.'])
        e_score_string = split_line[10]
        e_score = float(e_score_string)

        #Mark Gene as Present if Hit Exists over Threshold Value
        if e_score <= cutoff_float:
            if options['print_matches'] and not genes[gene]:
                analysis += 'Match: %s %s %s %s\n' % (sequence, BUSCO_sequence, gene, e_score_string)
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
    analysis += 'Percent BUSCO Genes Present: %s\n' % percent

    headers = ['Analys.', 'Cutoff', 'Expect.', 'Found', 'Missing', 'Total', 'Percent']
    data_grid = ['BUSCO', options['evalue_cutoff'], expected_gene_count, found_gene_count, 
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

            
    
    
    

