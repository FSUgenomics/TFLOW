#TFLOW Component: Dynamic Local Python Settings Based on File: "TFLOW_Local_Settings.dat"
#
#Dan Stribling
#Florida State University
#Center for Genomics and Personalized Medicine
#Version 0.9, 04/20/2015
#Project URL: http://www.github.com/fsugenomics/tflow

import os.path

ALLOWED_LOCAL_SETTINGS = ['TRIMMOMATIC_LOCATION',
                          'TRIMMOMATIC_EXEC',
                          'TRINITY_LOCATION',
                          'TRINITY_EXEC',
                          'CAP3_LOCATION',
                          'CAP3_EXEC',
                          'CEGMA_FILE',
                          'BUSCO_LOCATION',
                          'BLAST_LOCATION',
                          'BLAST_EXEC',
                          'MAKE_BLAST_DB_LOCATION',
                          'MAKE_BLAST_DB_EXEC',
                          ]

SETTINGS_DIRECTORY = os.path.dirname(os.path.realpath(__file__))
if os.path.isfile(os.path.join(SETTINGS_DIRECTORY, 'TFLOW_Local_Settings.dat')):
    LOCAL_SETTINGS_FILE = os.path.join(SETTINGS_DIRECTORY, 'TFLOW_Local_Settings.dat')

elif os.path.isfile(os.path.join(SETTINGS_DIRECTORY, '..', 'TFLOW_Local_Settings.dat')):
    LOCAL_SETTINGS_FILE = os.path.join(SETTINGS_DIRECTORY, '..', 
                                       'TFLOW_Local_Settings.dat')

elif os.path.isfile(os.path.join(SETTINGS_DIRECTORY, '..', '..' 
                               'TFLOW_Local_Settings.dat')):
    LOCAL_SETTINGS_FILE = os.path.join(SETTINGS_DIRECTORY, '..', '..', 
                                       'TFLOW_Local_Settings.dat')
else:
    LOCAL_SETTINGS_FILE = None

if LOCAL_SETTINGS_FILE:
    local_settings_temp = open(LOCAL_SETTINGS_FILE, 'r')
    for line in local_settings_temp:
        if not line.split() or line.startswith(('#', '!')):
            continue
        split_line = line.split()
        setting = split_line[0]
        value = split_line[1]
        if setting in ALLOWED_LOCAL_SETTINGS:
            exec("%s = '%s'" % (setting, value))

    local_settings_temp.close()
