#TFLOW Segment: Assemble RNA-Seq reads into transcripts with the Trinity Assembler
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import os.path
import sys
import subprocess
from collections import OrderedDict

if __name__ == "__main__" or __package__ is None:
    sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../'))
    import tflow.segments
    __package__ = "tflow.segments"

from .. import local_settings
from .parser_class import OutputParser
from ..fasta import check_N50_in_place
from ..util import (print_exit, print_error, print_warning, write_file, write_report, 
                    read_file_list, ensure_FASTQ_GZ, ensure_FASTA_GZ)

if hasattr(local_settings, 'TRINITY_LOCATION'):
    TRINITY_LOCATION = local_settings.TRINITY_LOCATION
else:
    TRINITY_LOCATION = ''

if hasattr(local_settings, 'TRINITY_EXEC'):
    TRINITY_EXEC = local_settings.TRINITY_EXEC
else:
    TRINITY_EXEC = os.path.join(TRINITY_LOCATION, 'Trinity')


JOB_TYPE = 'Trinity'
PROGRAM_URL = 'http://trinityrnaseq.github.io/'
SEGMENT_FOR_VERSION = 'trinityrnaseq_r20140717'
COMMAND = TRINITY_EXEC
COMMAND_LIST = [COMMAND]
TEST_COMMAND = '-help'
OUT_FILE = 'Trinity.out'
OUT_SEQUENCE_FILE = 'Trinity.fasta'
MILESTONES = ['Jellyfish',
              'Inchworm',
              'TIMING KMER_DB_BUILDING',
              'TIMING PRUNING',
              'TIMING CONTIG_BUILDING',
              'Chrysalis',
              'Chrysalis: GraphFromFasta',
              'Counting k-mers...',
              'done, assigning k-mers...',
              'Phase 1: Collecting candidate weldmers',
              '...done Phase 1',
              'Phase 2: Reclustering iworm contigs',
              'done bubbling',
              'Chrysalis: ReadsToTranscripts',
              'Chrysalis initial stage completed',
              'Chrysalis: QuantifyGraph',
              'Butterfly',
              'All commands completed',
              'Butterfly assemblies are written',
              ]
TERMINAL_FLAGS = ['Trinity Job Complete']
FAILURE_FLAGS = ['Exiting Early...',
                 'Traceback',
                 'Not Found']
DEFAULT_SETTINGS = {'is_paired_reads':True,
                    'read_type':'fq',
                    'max_memory':'10G',
                    'max_CPU':'4',
                    'output_dir':'Trinity_Assembly',
                    'out_sequence_file':OUT_SEQUENCE_FILE,
                    'min_contig_length':'200',

                    #Read Parsing Settings (Trimmomatic)
                    'left_read_indicator':'1P',
                    'right_read_indicator':'2P',
                    #Read Parsing Settings (Illumina)
                    #'left_read_indicator':'_R1',
                    #'right_read_indicator':'_R2',

                    #Optional Trinity Settings:
                    #'--genome':'genome.fa',
                    #'--jaccard_clip':'FLAG',
                    #'--SS_lib_type':[STRING],
                    #'--normalize_reads':'FLAG',
                    #'--full_cleanup':'FLAG',

                    #Optional Advanced Trinity Settings:
                    #'--prep':'FLAG',
                    #'--full_cleanup_ET':'FLAG',
                    #'--no_cleanup;:'FLAG',
                    #'--min_kmer_cov':[INT],
                    #'--inchworm_cpu':[INT],
                    #'--no_run_inchworm':'FLAG',
                    #'--max_reads_per_Graph':[INT],
                    #'--min_glue':[INT],
                    #'--no_run_chrysalis':'FLAG',
                    #'--no_run_quantifygraph':'FLAG',
                    #'--chrysalis_output':[STRING],
                    #'--no_bowtie':'FLAG',
                    #'--bfly_opts':[STRING],
                    #'--PasaFly':'FLAG',
                    #'--CuffFly':'FLAG',
                    #'--group_pairs_distance':[INT],
                    #'--path_reinforcement_distance':[INT],
                    #'--no_triplet_lock':'FLAG',
                    #'--extended_lock':'FLAG',
                    #'--NO_EM_REDUCE':'FLAG',
                    #'--no_path_merging':'FLAG',
                    #'--min_per_id_same_path':[INT],
                    #'--max_diffs_same_path':[INT],
                    #'--max_internal_gap_same_path':[INT],
                    #'--bflyHeapSpaceMax':[STRING],
                    #'--bflyHeapSpaceInit':[STRING],
                    #'--bflyGCThreads':[INT],
                    #'--bflyCPU':[INT],
                    #'--bflyCalculateCPU':'FLAG',
                    #'--no_run_butterfly':'FLAG',
                    #'--bfly_jar':[STRING],
                    #'--normalize_max_read_cov':[INT],
                    #'--normalize_by_read_set':'FLAG',
                    #'--genome_guided_max_intron':[INT],
                    #'--genome_guided_use_bam':[STRING],
                    #'--genome_guided_min_coverage':[INT],
                    #'--genome_guided_min_reads_per_partition':[INT],
                    #'--genome_guided_CPU':[INT],
                    #'--genome_guided_sort_buffer':[STRING],
                    #'--GMAP_CPU':[INT],
                    #'--genome_guided_just_prep':'FLAG',
                    #'--grid_conf_file':[STRING],
                     
                    #TFLOW Trinity Settings
                    'command':COMMAND,
                    'command_list':COMMAND_LIST,
                    'test_command':TEST_COMMAND,
                    'program_URL':PROGRAM_URL,
                    'segment_for_version':SEGMENT_FOR_VERSION,
                    #TFLOW Writing Defaults, Used if Global Not Set
                    'write_report':True,
                    'write_command':True,
                    }

REQUIRED_SETTINGS = ['command_list', 'is_paired_reads', 'read_type', 'max_memory', 'write_report',
                     'write_command']

REQUIRED_ANALYSIS_SETTINGS = ['working_directory', 'output_dir', 'out_sequence_file', 
                              'write_report']

ADVANCED_OPTIONS = ['--genome', '--jaccard_clip','--SS_lib_type','--normalize_reads',
                    '--full_cleanup', '--prep', '--full_cleanup_ET', '--no_cleanup',
                    '--min_kmer_cov', '--inchworm_cpu', '--no_run_inchworm', 
                    '--max_reads_per_Graph', '--min_glue', '--no_run_chrysalis', 
                    '--no_run_quantifygraph', '--chrysalis_output', '--no_bowtie', 
                    '--bfly_opts', '--PasaFly', '--CuffFly', '--group_pairs_distance', 
                    '--path_reinforcement_distance', '--no_triplet_lock', '--extended_lock',
                    '--NO_EM_REDUCE', '--no_path_merging', '--min_per_id_same_path', 
                    '--max_diffs_same_path', '--max_internal_gap_same_path', 
                    '--bflyHeapSpaceMax', '--bflyHeapSpaceInit', '--bflyGCThreads', '--bflyCPU'
                    '--bflyCalculateCPU', '--no_run_butterfly', '--bfly_jar', 
                    '--normalize_max_read_cov', '--normalize_by_read_set', 
                    '--genome_guided_max_intron', '--genome_guided_use_bam', 
                    '--genome_guided_min_coverage', '--genome_guided_min_reads_per_partition',
                    '--genome_guided_CPU', '--genome_guided_sort_buffer', '--GMAP_CPU',
                    '--genome_guided_just_prep', '--grid_conf_file' ]


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

    full_out_sequence_file = os.path.join(options['working_directory'], 
                                          options['output_dir'],
                                          options['out_sequence_file'])
    if not os.path.isfile(full_out_sequence_file):
        print_warning('Expected Trinity Output File: %s ' % full_out_sequence_file
                      + 'Cannot Be Found.' )
        return ''

    results = check_N50_in_place(full_out_sequence_file, fail_exit=False, 
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

def test(options, silent=False):
    try:
        output = subprocess.check_output(options['command_list'] + [options['test_command']])
    except subprocess.CalledProcessError as p_exception:
        if int(p_exception.returncode) == 255:
            output =  p_exception.output
            if silent:
                return True
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
            print_exit('Required Option: %s for Trinity not given.' % required_option)

    #Paired Input Reads
    if options['is_paired_reads']:
        #Attempt to guess pairing for grouped reads
        
        if any(x in options for x in ['all_reads', 'all_reads_list']):
            all_reads = []
            if 'all_reads' in options:
                raw_all_reads = options['all_reads']
            elif 'all_reads_list' in options:
                if not os.path.isfile(options['all_reads_list']):
                    print_exit('All Reads List File: %s' % options['all_reads_list']
                               + ' Not Found!')
                raw_all_reads = read_file_list(options['all_reads_list'])
           
            left_reads = []
            right_reads = []
            for read in raw_all_reads:
                full_read = os.path.join(options['working_directory'], read)
                all_reads.append(full_reads)
                if options['left_read_indicator'] in full_read:
                    left_reads.append(full_read)
                elif options['right_read_indicator'] in full_read:
                    right_reads.append(full_read)
                else:
                    print_exit('Pairing of Read %s Cannot Be Established' % full_read
                               + 'Using Pairing Indicators '
                               + 'Left: "%s" ' % options['left_read_indicator']
                               + ' and Right: %s.\n' % options['right_read_indicator']
                               + 'Please set indicators.')

        elif all(x in options for x in ['left_reads', 'right_reads']):
            left_reads = []
            right_reads = []

            if isinstance(options['left_reads'], str):
                raw_left_reads = [options['left_reads']]
            else:
                raw_left_reads = options['left_reads']
            if isinstance(options['right_reads'], str):
                raw_right_reads = [options['right_reads']]
            else:
                raw_right_reads = options['right_reads']
            
            for read in raw_left_reads:
                left_reads.append(os.path.join(options['working_directory'], read))
            for read in raw_right_reads:
                right_reads.append(os.path.join(options['working_directory'], read))
            all_reads = left_reads + right_reads

        elif all(x in options for x in ['left_reads_list', 'right_reads_list']):
            if not os.path.isfile(options['left_reads_list']):
                print_exit('All Reads List File: %s' % options['all_reads_list']
                           + ' Not Found!')
            if not os.path.isfile(options['right_reads_list']):
                print_exit('All Reads List File: %s' % options['all_reads_list']
                           + ' Not Found!')
            raw_left_reads = read_file_list(options['left_reads_list'])
            raw_right_reads = read_file_list(options['right_reads_list'])
            left_reads = []
            right_reads = []
            for read in raw_left_reads:
                left_reads.append(os.path.join(options['working_directory'], read))
            for read in raw_right_reads:
                right_reads.append(os.path.join(options['working_directory'], read))
            all_reads = left_reads + right_reads

        else:
            print_exit('Required Options: left_reads, right_reads, for input'
                       + ' of paired end reads not found.')

    #Unpaired Input Reads
    else:
        single_reads = []
        if 'single_reads' in options:       
            for read in single_reads:
                single_reads.append(os.path.join(options['working_directory'], read))
        elif 'single_reads_list' in options:
            if not os.path.isfile(options['single_reads_list']):
                print_exit('Single Reads List File: %s' % options['single_reads_list']
                           + ' Not Found!')
            for read in read_file_list(options['single_reads_list']):
                single_reads.append(os.path.join(options['working_directory'], read))
        
        all_reads = single_reads

    #Ensure reads exist and are of correct type.
    for read in all_reads:
        if options['read_type'] == 'fq':
            ensure_FASTQ_GZ(read)
        elif options['read_type'] == 'fa':
            ensure_FASTA_GZ(read)
        else:
            print_exit('Provided read_type value %s not "fq" or "fa"' % options['read_type'])           

    print ''
    if options['is_paired_reads']:
        if len(left_reads) != len(right_reads):
            print_exit('Number of Left Reads Does Not Equal Number of Right Reads.')

        print 'Using Paired End Reads:'
        print '  Left Reads:'
        for read in left_reads:
            print '  --' + read
        print ''
        print '  Right Reads:'
        for read in right_reads:
            print '  --' + read
    else:
        print 'Using Unpaired Reads:'
        for read in single_reads:
            print '  --' + read
    print ''
    

    
    command_list = options['command_list'][:]
    command_list += ['--seqType', options['read_type']]
    command_list += ['--JM', options['max_memory']]
    if options['is_paired_reads']:
        command_list += ['--left'] + left_reads
        command_list += ['--right'] + right_reads
    else:
        command_list += ['--single'] + single_reads
    command_list += ['--CPU', options['max_CPU']]
    if 'output_dir' in options:
        command_list += ['--output', options['output_dir']]
    if 'min_contig_length' in options:
        command_list += ['--min_contig_length', options['min_contig_length']]

    for advanced_option in ADVANCED_OPTIONS:
        if advanced_option in options:
            if options[advanced_option] in [None, 'FLAG']:
                command_list.append(advanced_option)
            else:
                command_list += [advanced_option, str(options[advanced_option])]

    command = ' '.join(command_list)

    if options['write_command']:
        command_file = os.path.join(options['project_directory'],
                                    options['job_type'] + '.auto.sh')
        write_file(command_file, '#!/bin/sh\n' +command)
    print 'Running Command:\n    ' + command

    sys.stdout.flush()

    try:
        process = subprocess.Popen(command_list, stdout=sys.stdout, stderr=sys.stderr,
                                   cwd=options['project_directory'])
        process.wait()
        sys.stdout.flush()

    except KeyboardInterrupt:
        if __name__ != '__main__' and options['is_pipe']:
            sys.stdout = terminal_out
            sys.stderr = terminal_error
            out_file_stream.close()

        print 'Killing Trinity Process.'
        process.kill()
        raise

    expected_output = os.path.join(options['working_directory'], options['output_dir'],
                                   'Trinity.fasta')
    if os.path.isfile(expected_output):
        print 'Expected Output Trinity File Found:'
        print ' ', expected_output
    else:
        print_error(['Expected Output Trinity File:', '  ' + expected_output, 'Not Found!'], 1)

    analyze(options)
    
    print ''
    print 'Trinity Job Complete'

    if __name__ != '__main__' and options['is_pipe']:
        sys.stdout = terminal_out
        sys.stderr = terminal_error
        out_file_stream.close()


