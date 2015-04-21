#TFLOW Segment_Pipe: De Novo Assembly of RNA-Seq Reads into Transcript Sequences.
#
#Steps:
#Make_Read_Lists: Parses reads into lists based on provided paramaters
#Trimmomatic:     Trim reads based on given quality settings
#Trinity:         Assemble reads into transcript sequences
#CAP3:            Further assemble output transcripts into longer sequences
#CEGMA_Analysis:  Analyze gene recapture of CEGMA core eukaryotic genes
#BUSCO_Analysis:  Analyze gene recapture of BUSCO benchmark genes
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

from collections import OrderedDict
import os.path

steps = OrderedDict()
TRINITY_DIR = 'Trinity_Assembly'
CAP3_DIR = 'CAP3'
steps['Make_Read_Lists'] = {'raw_left_reads_list':'raw_left_reads_list',
                            'raw_right_reads_list':'raw_right_reads_list',
                            'raw_single_reads_list':'raw_unpaired_reads_list',
                            }

steps['Trimmomatic'] = {'file_list_name':'trim_files',
                        'out_dir':'Trimmed_Data',
                        'raw_left_reads_list':'raw_left_reads_list',
                        'raw_right_reads_list':'raw_right_reads_list',
                        'raw_single_reads_list':'raw_unpaired_reads_list',
                        'left_reads_list':'trimmed_left_reads_list',
                        'right_reads_list':'trimmed_right_reads_list',
                        'single_reads_list':'trimmed_unpaired_reads_list',
                        }

steps['Trinity'] = {'output':TRINITY_DIR,
                    'left_reads_list':'trimmed_left_reads_list',
                    'right_reads_list':'trimmed_right_reads_list',
                    'single_reads_list':'trimmed_unpaired_reads_list',
                    }

steps['CAP3'] = {'working_directory':CAP3_DIR,
                 'relative_input_file':os.path.join(TRINITY_DIR, 'Trinity.fasta'),
                 }

steps['CEGMA_Analysis'] = {'rel_input_analysis_file':(os.path.join(CAP3_DIR, 
                                                                   'Trinity.fasta.cap.combined')),
                           'working_directory':'CEGMA_Analysis',
                           'copy_input':True,
                           }

steps['BUSCO_Analysis'] = {'BUSCO_type':'vertebrata',
                           'rel_input_analysis_file':(os.path.join(CAP3_DIR, 
                                                                   'Trinity.fasta.cap.combined')),
                           'working_directory':'BUSCO_Analysis',
                           'copy_input':True,
                           }
