#TFLOW Segment_Pipe: Assemble Sequences in Input File with CAP3 Assembler, Then Evaluate 
#Gene recapture
#
#Steps:
#CAP3:            Further assemble output transcripts into longer sequences
#Package:         Copy and Zip Final Sequence Output
#CEGMA_Analysis:  Analyze gene recapture of CEGMA core eukaryotic genes
#BUSCO_Analysis:  Analyze gene recapture of BUSCO benchmark genes
#Summary:         Create Summary Report of Results
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

from collections import OrderedDict

steps = OrderedDict()
steps['CAP3'] = {'working_directory':'',
                 'write_result_name':True,
                 }

steps['Package'] = {'result_name_file':'CAP3.auto.result_name',
                    }

steps['CEGMA_Analysis'] = {'working_directory':'CEGMA_Analysis',
                           'result_name_file':'CAP3.auto.result_name',
                           }
steps['BUSCO_Analysis'] = {'working_directory':'BUSCO_Analysis',
                           'result_name_file':'CAP3.auto.result_name',
                           }

steps['Summary'] = {}
