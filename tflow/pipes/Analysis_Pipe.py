#TFLOW Segment_Pipe: Statistical Analysis and Gene recapture analysis using BUSCO and 
#    CEGMA databases.
#
#Steps:
#Stat_Analysis:   Perform Statistical Analysis on Sequence File
#CEGMA_Analysis:  Analyze gene recapture of CEGMA core eukaryotic genes
#BUSCO_Analysis:  Analyze gene recapture of BUSCO benchmark genes
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

from collections import OrderedDict

steps = OrderedDict()

steps['Stat_Analysis'] = {'copy_input_file':False,
                          }

steps['CEGMA_Analysis'] = {'working_directory':'CEGMA_Analysis',
                           'copy_input_file':True,
                           }

steps['BUSCO_Analysis'] = {'working_directory':'BUSCO_Analysis',
                           'copy_input_file':True,
                           }
steps['Summary'] = {}
