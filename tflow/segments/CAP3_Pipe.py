#TFLOW Segment_Pipe: Assemble Sequences in Input File with CAP3 Assembler, Then Evaluate 
#Gene recapture
#
#Steps:
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
steps['CAP3'] = {'working_directory':''
                 }

steps['CEGMA_Analysis'] = {'working_directory':'CEGMA_Analysis',
                           }

steps['BUSCO_Analysis'] = {'BUSCO_type':'vertebrata',
                           'working_directory':'BUSCO_Analysis',
                           }
