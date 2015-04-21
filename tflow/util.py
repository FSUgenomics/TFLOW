#TFLOW Component: Commonly Used Printing, Formatting, and Other Functions
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import sys
import os.path
import subprocess

# --- Global Constants ---

BOOL = [True, False]
RIGID_BOOL = {'F':False, 'FALSE':False, 'False':False, 'T':True, 'TRUE':True, 'True':True}
FLEXIBLE_FALSE = ['False', 'false', 'F', 'f', '0']
FLEXIBLE_TRUE = ['True', 'true', 'T', 't', '1']
FLEXIBLE_BOOL = FLEXIBLE_FALSE + FLEXIBLE_TRUE
ACTION_NAMES = {'run':'Running', 'track':'Tracking', 'analyze':'Analyzing', 
                'read':'Reading', 'test':'Testing', 'track-analyze':'Tracking and Analyzing'}
DEFAULT_SETTINGS = {'TFLOW_Version':'0.9',
                    'write_analysis':True,
                    'write_settings':True,
                    'write_command':True,
                    'write_times':True,
                    'is_pipe':False,
                    'verbose':False,
                    'overwrite':False,
                    'print_test_output':False,
                    }

# --- Output Functions ---
def print_except(message):
    print >> sys.stderr, message
    raise Exception(message)

def print_exit(message, code=1):
    if isinstance(message, list):
        for piece in message:
            print piece
    else:
        print message

    if code == 0:
        print 'Exiting...'
    else:
        print 'Exiting Early...'
    print ''
    sys.exit(code)


def print_error(message, exit_code=-1):
    if isinstance(message, list):
        for piece in message:
            print >> sys.stderr, 'ERROR:', piece
    else:
        print >> sys.stderr, 'ERROR:', message
    sys.stderr.flush()
    if exit_code > -1:
        sys.exit(exit_code)

def print_warning(message):
    if isinstance(message, list):
        for piece in message:
            print >> sys.stderr, 'WARNING:', piece
    else:
         print >> sys.stderr, 'WARNING:', message
    sys.stderr.flush()

def print_multi(*args):
    for arg in args:
        print arg
    sys.stdout.flush()

def write_file(name, contents):
    f = open(name, 'w')
    f.write(contents)
    f.close()

def write_file_list(name, input_list):
    list_file = open(name, 'w')
    for item in input_list:
        list_file.write(str(item)+'\n')
    list_file.close()

def write_date_time(name, start=None):
    from datetime import datetime
    now_time = datetime.today()
    if start:
        contents  = 'Start Time: ' + str(start) + '\n'
        contents += 'End Time:   ' + str(now_time) + '\n'
        contents += 'Total Time: ' + str(now_time - start)+'\n'
        write_file(name, contents)
    else:
        write_file(name, 'Start Time: ' + str(now_time))

    return now_time




# --- Reading Functions ---

# - Read a List From a File
def read_file_list(file_name):
    return_list = []
    if not os.path.isfile(file_name):
        print_exit('List File %s Not Found!' % file_name)
    list_file = open(file_name, 'r')
    for line in list_file:
        return_list.append(line.strip())
    list_file.close()
    return return_list


# - Read Options From A File 
FULL_OPTIONS_FILES = ['options.dat', 'job_options.dat', 'project_options.dat']
def get_file_settings():
    file_options = {}

    if not any(os.path.isfile(file_name) for file_name in (FULL_OPTIONS_FILES + ['jobtype.dat'])):
        return file_options

    if os.path.isfile('jobtype.dat'):
        type_file = open('jobtype.dat', 'r')
        file_options['job_type'] = type_file.read().strip()
        type_file.close()

    #Find Options File Name
    for file_name in FULL_OPTIONS_FILES:
        if os.path.isfile(file_name):
            break

    if os.path.isfile(file_name):
        options_file = open('options.dat', 'r')
        for line in options_file:
            split_line = line.split()
            if not split_line or split_line[0].startswith(('#', '!')):
                continue

            key = split_line[0]
            if len(split_line) > 2:
                data = []
                for datum in split_line[1:]:
                    if datum.startswith(('#', '!')):
                        break
                    data.append(filter_boolean(datum))
                if len(data) == 1:
                    data = data[0]
            else:
                data = filter_boolean(split_line[1])

            if '.' in key:
                split_key = key.split('.')
                if len(split_key) > 2:
                    print_except('Problem!!! Settings Key: %s ' % key
                                 + 'Has More Than Two Variables!')
                job, key = split_key[0], split_key[1]
            else:
                job = None

            if job:
                if job not in file_options:
                    file_options[job] = {}
                file_options[job][key] = data
            else:
                file_options[key] = data

        options_file.close()
    return file_options


def read_process_output(process):
    for line in process.stdout:
        print line.rstrip()

# --- Formatting Functions ---

def lowercase(string):
    return (str(string).lower())

def flexible_boolean_string(in_var):
    if in_var not in FLEXIBLE_BOOL:
        return in_var
    if in_var in ['False', 'false', 'F', 'f', '0']:
        return False
    return True

def flexible_boolean(in_var):
    if not in_var or in_var in ['False', 'false', 'F', 'f', '0']:
        return False
    return True

def filter_boolean(in_var):
    if in_var in RIGID_BOOL:
        return RIGID_BOOL[in_var]
    return in_var

def SI_prefix(number):
    if number > 1000000000000:
        formatted_number = number/1000000000000
        prefix = 'T'
    elif number > 1000000000:
        formatted_number = number/1000000000
        prefix = 'G'
    elif number > 1000000:
        formatted_number = number/1000000
        prefix = 'M'
    elif number > 1000:
        formatted_number = number/1000
        prefix = 'k'
    else:
        formatted_number = number
        prefix = ''

    return (formatted_number, prefix)

def percent_string(numerator, denominator):
    if not denominator:
        return 'N/A%'
    return '{0:.3g}'.format((float(numerator)/float(denominator))*100)+'%'


def return_settings(options, message=None):
    return_string = ''
    if message:
        return_string += str(message) + '\n'
        print 'TFLOW Conductor Options:'

    max_len = len(max(options.keys(), key=len))
    for option in sorted(options.keys()):
        if isinstance(options[option], dict):
            return_string += '-- %s OPTION SUBSET:\n' % option
            sub_options = options[option]
            for sub_option in sorted(sub_options.keys()):
                max_sublen = len(max(sub_options.keys(), key=len))
                return_string += ('---- ' + sub_option.ljust(max_sublen) + ' '
                                  + str(sub_options[sub_option]) + '\n')

        else:
            return_string += '-- ' + option.ljust(max_len) + ' ' + str(options[option]) + '\n'

    return return_string

def print_settings(options, message=None):
    print return_settings(options, message=None)
    print ''

def write_settings(options, out_file, message=None):
    write_file(out_file, return_settings(options, message=None))

# --- Sequence File Utilities ---
def count_FASTA(file_name):
    return int((subprocess.check_output('grep -c "^>" %s' % file_name, 
                                        shell=True)).strip())

def count_FASTA_GZ(file_name):
    return int((subprocess.check_output('zcat %s | grep -c "^>"' % file_name, 
                                        shell=True)).strip())

def count_FASTA_all(file_name):
    if is_FASTA(file_name):
        return count_FASTA(file_name)
    elif is_FASTA_GZ(file_name):
        return count_FASTA_GZ(file_name)
    return 0

def count_FASTQ(file_name):
    return (int((subprocess.check_output('wc -l %s' % file_name, 
                                         shell=True)).split()[0]) / 4)

def count_FASTQ_GZ(file_name):
    return (int((subprocess.check_output('zcat %s | wc -l' % file_name, 
                                         shell=True)).split()[0]) / 4)

def count_FASTQ_all(file_name):
    if is_FASTQ(file_name):
        return count_FASTQ(file_name)
    elif is_FASTQ_GZ(file_name):
        return count_FASTQ_GZ(file_name)
    return 0

def is_FASTA(file_name):
    return file_name.endswith(('.fasta', '.fa', '.fas', '.fna', 'ffn', 'faa', 'frn'))

def is_FASTA_GZ(file_name):
    return file_name.endswith(('.fasta.gz', '.fa.gz', '.fas.gz', '.fna.gz', 'ffn.gz', 
                               'faa.gz', 'frn.gz')) 

def is_FASTQ(file_name):
    return file_name.endswith(('.fastq', '.fq'))

def is_FASTQ_GZ(file_name):
    return file_name.endswith(('.fastq.gz', '.fq.gz'))

def is_sequence_file(file_name):
    return (is_FASTA(file_name) or is_FASTA_GZ(file_name) or 
            is_FASTQ(file_name) or is_FASTQ_GZ(file_name))

def ensure_FASTA(file_name):
    if not os.path.isfile(file_name):
        print_exit('File: %s cannot be found.' % file_name)
    if not is_FASTA(file_name):
        print_exit('File: %s is not a *.fasta file.' % file_name)
    return True

def ensure_FASTA_GZ(file_name):
    if not os.path.isfile(file_name):
        print_exit('File: %s cannot be found.' % file_name)
    if not (is_FASTA(file_name) or is_FASTA_GZ(file_name)):
        print_exit('File: %s is not a *.fasta or *.fasta.gz file.' % file_name)
    return True

def ensure_FASTA_only_GZ(file_name):
    if not os.path.isfile(file_name):
        print_exit('File: %s cannot be found.' % file_name)
    if not is_FASTA_GZ(file_name):
        print_exit('File: %s is not a *.fasta.gz file.' % file_name)
    return True

def ensure_FASTQ(file_name):
    if not os.path.isfile(file_name):
        print_exit('File: %s cannot be found.' % file_name)
    if not is_FASTQ(file_name):
        print_exit('File: %s is not a *.fastq file.' % file_name)
    return True

def ensure_FASTQ_GZ(file_name):
    if not os.path.isfile(file_name):
        print_exit('File: %s cannot be found.' % file_name)
    if not (is_FASTQ(file_name) or is_FASTQ_GZ(file_name)):
        print_exit('File: %s is not a *.fastq or *.fastq.gz file.' % file_name)
    return True

def ensure_FASTQ_only_GZ(file_name):
    if not os.path.isfile(file_name):
        print_exit('File: %s cannot be found.' % file_name)
    if not is_FASTQ_GZ(file_name):
        print_exit('File: %s is not a *.fastq.gz file.' % file_name)
    return True


def count_sequences(contents, print_results=True, column_width=50):
    counts = []
    total_count = 0
    for file_name in contents:
        if is_FASTA(file_name):
            count = count_FASTA(file_name)

        elif is_FASTA_GZ(file_name):
            count = count_FASTA_GZ(file_name)

        elif is_FASTQ(file_name):
            count = count_FASTQ(file_name)

        elif is_FASTQ_GZ(file_name):
            count = count_FASTQ_GZ(file_name)

        else:
            print >> sys.stderr, ('File Format for File: %s Not Identified, ' % file_name
                                  + 'Assuming FASTA' )
            count = count_FASTA(file_name)

        if count != None:
            if print_results:
                print file_name.ljust(column_width), count
            counts.append((file_name, count))
            total_count += count

    if print_results and len(counts) > 1:
        print 'Total'.ljust(column_width), total_count
    counts.append(('Total', total_count))

    return counts



# --- Type Conversion Utilities ---
def string_to_boolean(string):
    if not string or string.lower in ['f', 'false']:
        return False
    else:
        return True

def ensure_list(in_variable):
    if isinstance(in_variable, list):
        return in_variable
    else:
        return [ in_variable ]





