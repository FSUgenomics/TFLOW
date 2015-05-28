#TFLOW Component: FASTA Sequence File Utilities and Database Class
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import os.path
import sys
import re

from collections import OrderedDict
from .util import print_except, print_exit, SI_prefix, percent_string, is_FASTA, is_FASTQ


#Fasta Database Class
class FASTA_DB():
    def __init__(self):
        self.reset_data()
        self.reset_settings()

    def reset_data(self):
        self.in_file = None
        self.out_file = None
        self.base_name = None
        self.sequence_type = None
        self.sequences = OrderedDict()
        self.header = ''

    def reset_settings(self):
        self.except_identical_name = False
        self.except_empty_contents = True
        self.except_file_exists = False
        self.except_extra_lines = True
        self.except_duplicate_match = True

    def set_lenient(self):
        self.except_identical_name = False
        self.except_empty_contents = False
        self.except_file_exists = False
        self.except_extra_lines = False
        self.except_duplicate_match = False

    def set_strict(self):
        self.except_identical_name = True
        self.except_empty_contents = True
        self.except_file_exists = True
        self.except_extra_lines = True
        self.except_duplicate_match = True

    def add_sequence(self, name, contents, line_number=-1):
        if name in self.sequences:
            if self.except_identical_name:
                print_except('ERROR: Repeated Sequence Header: %s' % name)
            else:
                for header in ['Repeat' + str(x) + '_' for x in range(1,100)]:
                    new_name = header +  name
                    if new_name not in self.sequences:
                        break
                print >> sys.stderr, ('ERROR: Repeated Sequence Header: '
                                      '%s   Added as %s' % (name, new_name))
                name = new_name                

        if not contents:
            if self.except_empty_contents:
                print_except('ERROR: Sequence %s ' % name
                             + 'from line: %s ' % str(line_number)
                             + 'has no detectable contents') 
            else:
                print >> sys.stderr, ('ERROR: Sequence %s ' % name
                                     + 'has no detectable contents')  

        self.sequences[name] = contents
    
    def read_file(self, file_name, reset=True):
        if not os.path.isfile(file_name):
            print 'File %s Not Found, No Changes Made.' % file_name
            return False

        if reset:
            self.reset_data()
            start_sequence_count = 0
        else:
            start_sequence_count = len(self.sequences)

        print 'Reading File: %s' % file_name
        f = open(file_name, 'rb')
        
        f_enumerator = enumerate(f, start=1)        
        
        try:
            (line_count, line) = f_enumerator.next()          
            while line.startswith(';'):
                self.header += line.lstrip(';')
                (line_count, line) = f_enumerator.next()

            if self.header:
                print 'Sequence File Header:'
                print self.header.rstrip()

            if not line or not line.startswith('>'):
                print_except('FASTA File %s Formatted Incorrectly, ' % file_name
                             + 'Header Space Found at line %i.' % line_count)
            name = line.lstrip('>').strip()
            sequence_line = line_count
            (line_count, line) = f_enumerator.next()

        except StopIteration:
            print 'Empty Sequence File.'
            return 0

        if line.startswith('>'):
            print_except('FASTA File %s Formatted Incorrectly, First ' % file_name
                         + '"Contents" line at line %i is also a sequence header.' % line_count)

        contents = line.strip()

        for (line_count, line) in f_enumerator:
            line = line.strip()
            if not line:
                if self.except_extra_lines:
                    print_except('FASTA File %s Formatted Incorrectly, Blank Line ' % file_name
                                 + 'Found at Line %i.' % line_count)
                else:
                    print >> sys.stderr, ('ERROR: Blank Line Found '
                                          + 'at Line %i.' % line_count)
                    continue                 

            if line.startswith('>'):
                #Add Previous Sequence
                self.add_sequence(name, contents, sequence_line)

                #Start New Sequence 
                name = line.lstrip('>').strip()
                sequence_line = line_count
                contents = ''
            else:
                contents += line.strip()

        #Add Last Sequence
        self.add_sequence(name, contents, sequence_line)

        if reset:
            self.in_file = file_name
            self.set_base_name(file_name)
         
        f.close()

        sequences_added = (len(self.sequences) - start_sequence_count)

        print '%i Sequences Read from File %s' % (sequences_added, file_name)
        return sequences_added


    def add_file(self, file_name):
        return read_file(self, file_name, reset=False)


    def print_sequences(self, space=True):
        if self.header:
            print self.header
            if space:
                print ''

        for name in self.sequences:
            print name
            print self.sequences[name]
            if space:
                print ''

    def get_details(self):
        shortest = (10000000000, 'None')
        longest = (0, 'None')
        for name in self.sequences.keys():
            sequence_length = len(self.sequences[name])
            if sequence_length < shortest[0]:
                shortest = (sequence_length, name)
            if sequence_length > longest[0]:
                longest = (sequence_length, name)

        report =  'Shortest Sequence: %i bp,\t%s\n' % shortest 
        report += 'Longest Sequence:  %i bp,\t%s\n' % longest

        return report

    def _OLD_Read_FromFASTQ(self, file_name, verify=True):
        if os.path.isfile(file_name):
            print 'Reading FASTQ File: %s' % file_name
            in_file = open(file_name, 'r')
            in_fileLines = in_file.readlines()
            in_file.close()

            self.reads=[]
            for i in range(0, len(in_fileLines), 4):
                readName = in_fileLines[i].strip('@\n')
                data = in_fileLines[i+1].strip()
                self.reads.append([readName, data])

            self.in_fileN = file_name
            self.set_base_name(file_name)

            print len(self.reads), 'Reads Read from File.'
            if verify:
                self.verify_reads()


        else:
            print 'File Not Found, No Changes Made.'
            return

    def set_base_name(self, file_name):
        if is_FASTA(file_name):
            self.base_name = (file_name.rstrip('fstar')).rstrip('.')
        elif is_FASTQ(file_name):
            self.base_name = (file_name.rstrip('fastq')).rstrip('.')
        else:
            self.base_name = file_name
        return self.base_name


    def write_file(self, file_name, overwrite=False):
        if os.path.exists(file_name) and not overwrite:
            message = 'Attempting to Write to File: %s, but it Already Exists!' % file_name
            if self.except_file_exists:
                read_except(message)
            else:
                print >> sys.stderr, 'ERROR: '+message+', No Changes Made.'
                return False

        print 'Writing Sequences to File: %s' % file_name
        out_file = open(file_name, 'w')
        if self.header:
            for line in self.header.splitlines(False):
                out_file.write(';' + line + '\n')

        for sequence_name in self.sequences:
            out_file.write('>'+sequence_name+'\n')
            formatted_sequence = '\n'.join(re.findall(r'.{1,60}', 
                                           self.sequences[sequence_name], 
                                           re.DOTALL))
            out_file.write(formatted_sequence+'\n')

        out_file.close()

        self.out_file = file_name
        print len(self.sequences), 'Sequences Written to file %s.' % file_name
        return True

    def N50(self):
        print 'Finding N50 Length of Sequences',
        print ''

        print 'Counting Sequence Lengths'

        sequence_length_counts = {} 
        counted_sequences = 0
        total_length = 0 

        for name in self.sequences:
            length = len(self.sequences[name])
            total_length += length
            counted_sequences += 1
            if length in sequence_length_counts:
                sequence_length_counts[length] += 1
            else:
                sequence_length_counts[length] = 1

        if counted_sequences != len(self.sequences):
            raise_except('PROBLEM!!! Counted_sequences != Num Sequences')

        print 'Finding Average Length...'
        average_length = total_length/counted_sequences
        sequence_lengths = sequence_length_counts.keys()
        sequence_lengths.sort()

        (formatted_number, prefix) = SI_prefix(total_length)
        formatted_total_length = str(formatted_number) + ' ' + prefix + 'bp'

        print 'Finding Median...'
        median_index = counted_sequences/2
        encountered = 0
        for length in sequence_lengths:
            encountered += sequence_length_counts[length]
            if encountered >= median_index:
                median_length = length
                break

        print 'Finding N50...'
        weighted_counts = {}
        total_weighted_counts = 0

        for length in sequence_lengths:
            weighted_count =  length * sequence_length_counts[length]
            weighted_counts[length] = weighted_count
            total_weighted_counts += weighted_count

        weighted_half = total_weighted_counts/2
        encountered = 0
        for length in sequence_lengths:
            encountered += weighted_counts[length]
            if encountered >= weighted_half:
                n50_length = length
                break

        print ''
        results = 'Results:\n'
        results += 'Number of Sequences: %i\n' % counted_sequences
        results += 'Total Sequence Length: %s\n' % formatted_total_length 
        results += 'Average Sequence Length: %i\n' % average_length
        results += 'Range: %i to %i\n' % (sequence_lengths[0], sequence_lengths[-1])
        results += 'Median Length: %i\n' % median_length
        results += 'N50 Length: %i\n' % n50_length
        results += 'Count, Len, Av.Len, SRange, ERange, Median, N50\n'
        results += '%i, %s, %i, %i, %i, %i, %i\n' % (counted_sequences, formatted_total_length, 
                                                     average_length, sequence_lengths[0], 
                                                     sequence_lengths[-1], median_length, n50_length)
        print results
        return results   
   

    def Cull(self, matches, ensure_match=True):
        #Prefered matching input is dictionary, but an input list 
        #will be converted.
        print 'Culling Sequences...'
 
        starting_length = len(self.sequences)
        print 'Starting Sequences:', starting_length
        print 'Match Terms:', len(match_list)

        if isinstance(matches, list):
            print 'Creating Match Dictionary:'
            match_dict = OrderedDict()
            for match in matches:
                if match not in match_dict:
                    match_dict[match] = None
                else:
                    if self.except_duplicate_match:
                        print_except('Match Item: %s is a duplicate!')
                    else:
                        print >> sys.stderr, 'ERROR: Match Item: %s is a duplicate!' % match

            matches = match_dict

        culled_sequences = OrderedDict()

        for name in self.sequences:
            if name not in matches:
                culled_sequences[name] = self.sequences[name]
                del(self.sequences[name])

        if ensure_match:
            for match in matches:
                if match not in self.sequences:
                    print_except('Match ID: %s Not found in database!' % match)

        final_length = len(self.sequences)

        if (final_length + len(self.sequences)) != starting_length:
            print_except('PROBLEM: Culled + Remaining != Starting!')

        print 'Sequence Culling Complete.'
        print 'Culled Sequences', len(culled_sequences)
        print 'Final Sequences:', final_length
        print percent_string(final_length, starting_length), ' of Sequences Retained.'


    def label_prefix(self, prefix, connector='_'):
        start_length = len(self.sequences)
        for name in self.sequences.keys():
            self.sequences[prefix+connector+name]=self.sequences[name]
            del(self.sequences[name])

        if start_length != len(self.sequences):
            print_except('PROBLEM: Num Labeled Sequences != Num Initial Sequences.')  


    def _OLD_SpecialLabel(self, tissue, destinations):
        if not self._HasSearchDB():
            print ''
            self.MakeSearchDB()
            print ''

        for read in self.reads:
            name = read[0].split()[0]
            if name not in destinations:
                print 'Sequence %s Not Found in Desitnation Dictionary' % name
                sys.exit(1)
            
            if destinations[name].rstrip('0123456789') == 'Contig':
                destination = destinations[name]
            else:
                destination = 'NA'
            read[0] = tissue+'|'+destination+'|'+read[0]
                              

    def trim_names(self, length=1000):
        new_sequences = OrderedDict()

        start_length = len(self.sequences)
        for name in self.sequences:           
            if len(name) >= length:
                new_sequences[name[:length]] = self.sequences[name]
            else:
                new_sequences[name] = self.sequences[name]
            del(self.sequences[name])

        if start_length != len(new_sequences):
            print_except('PROBLEM: Num Labeled Sequences != Num Initial Sequences.')  

        self.sequences = new_sequences          


    # --Python Magic Methods--
    # Enable iteration over class instance
    def __iter__(self):
        return self.sequences.keys().__iter__()

    # Equality comparison, based on whether Annotation instances have same name.
    def __eq__(self, other):
        if not isistance(other, FASTA_DB):
            return False
        else:
            return self.sequences == other.sequences

    # Return String Describing instance
    def __unicode__(self):
        if self.basename:
            name = self.basename
        else:
            name = 'EMPTY'
        return '<FASTA_DB Instance:%s>' % name

    # Length datum, return number of sequences in db
    def __len__(self):
        return len(self.sequences)

    # Boolean comparison, returns True if instance contains records not empty.
    def __nonzero__(self):
        return bool(self.sequences)

    # Item in container:
    def __contains__(self, name):
        return name in self.sequences

    # Retrieve annotation by name.
    def __getitem__(self, name):
        return self.sequences[name]


#Checks a FASTA file for formatting issues.
def check_FASTA(file_name):

    if not os.path.isfile(file_name):
        print 'FASTA File: %s Does Not Exist.' % file_name
        return False

    database = FASTA_DB()
    database.read_file(file_name)
    return True

def check_N50(file_name):
    if not os.path.isfile(file_name):
        print 'File %s Does Not Exist.' % file_name
        return False

    sequences = FASTA_DB()
    sequences.read_file(file_name)
    print ''
    sequences.N50()
    return True


def check_N50_in_place(file_name, fail_exit=True, return_report=False, return_report_dict=False):
    analysis = ''
    if not os.path.isfile(file_name):
        print 'File %s Not Found, No Changes Made.' % file_name
        return analysis

    sequence_length_counts = {} 
    counted_sequences = 0

    print 'Reading File: %s' % file_name
    f = open(file_name, 'r')
    
    f_enumerator = enumerate(f, start=1)        

    (line_count, line) = f_enumerator.next()
    while line.startswith(';'):
        (line_count, line) = f_enumerator.next()

    if not line or not line.startswith('>'):
        if fail_exit:
            print_exit('FASTA File %s Formatted Incorrectly, ' % file_name
                         + 'Header Space Found at line %i.' % line_count)
        else:
            analysis += ('FASTA File %s Formatted Incorrectly, ' % file_name
                         + 'Header Space Found at line %i.' % line_count)
            return analysis
        
    (line_count, line) = f_enumerator.next()
    if line.startswith('>'):
        if fail_exit:
            print_exit('FASTA File %s Formatted Incorrectly, First "Contents" ' % file_name
                         + 'line at line %i is also a sequence header.' % line_count)
        else:
            analysis += ('FASTA File %s Formatted Incorrectly, First "Contents" ' % file_name
                         + 'line at line %i is also a sequence header.' % line_count)
            return analysis

    content_length = len(line.strip())

    for (line_count, line) in f_enumerator:
        line = line.strip()
        if not line:
            if fail_exit:
                print_exit('FASTA File %s Formatted Incorrectly, ' % file_name
                             + 'Blank Line Found at Line %i.' % line_count)
            else:
                analysis += ('FASTA File %s Formatted Incorrectly, ' % file_name
                             + 'Blank Line Found at Line %i.' % line_count)

        if line.startswith('>'):
            #Add Previous Sequence Length
            if content_length in sequence_length_counts:
                sequence_length_counts[content_length] += 1
            else:
                sequence_length_counts[content_length] = 1
            counted_sequences += 1

            #Start New Sequence 
            content_length = 0

        else:
            content_length += len(line.strip())

    #Add Last Sequence
    if content_length in sequence_length_counts:
        sequence_length_counts[content_length] += 1
    else:
        sequence_length_counts[content_length] = 1
    counted_sequences += 1

    f.close()

    print '%i Sequence Lengths Read from File %s' % (counted_sequences, file_name)

    print 'Finding Total and Average Length...'
    total_length = 0 
    for length in sequence_length_counts:
        total_length += length * sequence_length_counts[length]

    average_length = total_length/counted_sequences
    sequence_lengths = sequence_length_counts.keys()
    sequence_lengths.sort()

    (formatted_number, prefix) = SI_prefix(total_length)
    formatted_total_length = str(formatted_number) + ' ' + prefix + 'bp'

    print 'Finding Median...'
    median_index = counted_sequences/2
    encountered = 0
    for length in sequence_lengths:
        encountered += sequence_length_counts[length]
        if encountered >= median_index:
            median_length = length
            break


    print 'Finding N50...'
    weighted_counts = {}
    total_weighted_counts = 0

    for length in sequence_lengths:
        weighted_count =  length * sequence_length_counts[length]
        weighted_counts[length] = weighted_count
        total_weighted_counts += weighted_count

    weighted_half = total_weighted_counts/2
    encountered = 0
    for length in sequence_lengths:
        encountered += weighted_counts[length]
        if encountered >= weighted_half:
            n50_length = length
            break

    print ''
    analysis += 'Results:\n'
    analysis += 'Number of Sequences: %i\n' % counted_sequences
    analysis += 'Total Sequence Length: %s\n' % formatted_total_length 
    analysis += 'Average Sequence Length: %i\n' % average_length
    analysis += 'Range: %i to %i\n' % (sequence_lengths[0], sequence_lengths[-1])
    analysis += 'Median Length: %i\n' % median_length
    analysis += 'N50 Length: %i\n' % n50_length
    analysis += '\n'
    analysis += 'Tab-Separated Results:\n'
    

    keys = ['Count', 'Len', 'Av.Len', 'SRange', 'ERange', 'Median', 'N50']
    values = [str(x) for x in [counted_sequences, formatted_total_length, average_length, 
                               sequence_lengths[0], sequence_lengths[-1], median_length, 
                               n50_length]]


    report_dict = dict(zip(keys, values))
    report_dict['report_type'] = 'sequence'
    report = '\t'.join(keys) + '\n'
    report += '\t'.join(values) + '\n'

    analysis += report
                        
    print analysis
    if return_report or return_report_dict:
        if return_report_dict:
            return analysis, report_dict
        else:
            return analysis, report
    else:
        return analysis


def label_sequences(in_file_name, label, out_file_name=None):
    if not os.path.isfile(in_file_name):
        print 'File %s Does Not Exist.' % in_file_name
        return False

    print 'Adding Label: %s to Sequences in File: %s' % (label, in_file_name)
    print ''
    database = FASTA_DB()
    database.read_file(in_file_name)
    print ''
    database.label_prefix(label)
    if not out_file_name:
        out_file_name = in_file_name + '.labeled'

    database.write_file(out_file_name)
    return True


def details(in_file_name):
    if not os.path.isfile(in_file_name):
        print 'File %s Does Not Exist.' % in_file_name
        return False

    database = FASTA_DB()
    database.read_file(in_file_name)
    print ''
    details = database.get_details()
    print details
    return True
