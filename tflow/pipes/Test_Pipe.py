#TFLOW Test_Pipe: Pipe to Test All Supported TFLOW Segments
#
#Steps:
#Make_Read_Lists, Trimmomatic, Trinity, CAP3, CEGMA_Analysis, BUSCO_Analysis
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

from collections import OrderedDict
import os.path

steps = OrderedDict()
steps['Make_Read_Lists'] = {}
steps['Trimmomatic'] = {}
steps['Trinity'] = {}
steps['CAP3'] = {}
steps['Stat_Analysis'] = {}
steps['CEGMA_Analysis'] = {}
steps['BUSCO_Analysis'] = {}
steps['Package'] = {}
steps['Summary'] = {}
