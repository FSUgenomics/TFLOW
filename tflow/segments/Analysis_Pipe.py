#TFLOW Segment_Pipe: Gene recapture analysis using BUSCO and CEGMA databases.
#
#Steps:
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

steps['CEGMA_Analysis'] = {'working_directory':'CEGMA_Analysis',
                           'copy_input':True,
                           }

steps['BUSCO_Analysis'] = {'BUSCO_type':'vertebrata',
                           'working_directory':'BUSCO_Analysis',
                           'copy_input':True,
                           }
