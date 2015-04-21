

Welcome to TFLOW, a command-line pipeline for de novo transcriptome assembly and analysis.

 --- Get Me Started ---
Optional: Add TFLOW dir to your path with: "export PATH=[TFLOW_DIR]:$PATH"
Optional: Add TFLOW dir to your python path with: "export PYTHONPATH=[TFLOW_DIR]:$PYTHONPATH"
Run: "tflow.sh test -t Trinity_Pipe"

Then, for each component that can't be found:
    Install (if not already)
    in TFLOW_Local_Settings.dat:
        set [UTILITY]_LOCATION variable with executable conatining folder
    	OR: set [UTILITY]_EXEC variable of full executable path

When every component (that you care about) reports "Found!":
    Create a project directory
    Copy Desired Pipe Options File Template to Directory:
    	Ex: Copy options.dat.Trinity_Pipe_Template to My_Project/options.dat)
    Set Mandatory and Suggested Settings in Options File:
    	Ex: is_paired_reads <True/False>
	    label	    <tissue_name>
	    BUSCO_type	    <vertebrata>
    Run:
        Ex: "TFLOW_DIR/tflow.sh run --raw_reads MyData/*.gz 
	OR: "TFLOW_DIR/tflow.sh run --raw_left_reads MyData/*left_reads*.gz --raw_right_reads MyData/*right_reads*.gz 

    Then, See what's going on by running:
        "TFLOW_DIR/tflow.sh run -h"
        OR: "TFLOW_DIR/tflow.sh run --help"


 --- Operation ---
TFLOW 0.9 has 6 modes:
    run	      	    Run Job/Pipe
    track	    Track Job/Pipe
    analyze	    Analyze Job/Pipe
    read 	    Read Output of Job/Pipe
    test	    Test Job/Pipe
    print_settings  Print Current TFLOW Settings


 --- Design ---
The TFLOW pipeline is designed with several key features in mind
    Easy accessability of assembly utilities
    A smooth "Flow" of input/output data from one tool to the next
    Automated analysis of results
    Modular addition of new tools

 --- Options ---
The priority hierarchy of options is given by:
Comand Line > Options File > Global Defaults (Should Not Conflict) > Pipe Settings File 
     > Segment Defaults  


--- Version History ---
2015-04-20, V.0.9, Initial Beta Release
