[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3817475.svg)](https://doi.org/10.5281/zenodo.3817475)

Welcome to TFLOW, a command-line pipeline for de novo transcriptome assembly and analysis.

 --- Design ---
 
The TFLOW pipeline is designed with several key features in mind:
    Easy accessability of assembly utilities
    A smooth "Flow" of input/output data from one tool to the next
    Automated analysis of results
    Modular addition of new tools
    TFLOW is built on the Python2.7 language and thus requires Python2.7 as a dependency.
    Individual external components are required depending on the pipe utilized.


 --- Get Me Started ---

  - Setup
    1. Download and Unpack TFLOW to a directory of your choice.
        Ex: unzip ~/v0.9.2.zip (creates ~/TFLOW-0.9.1)
        
    2. Run test_TFLOW_setup.sh
        Ex: ~/TFLOW-0.9.1/test_TFLOW_setup.sh
            (equivalent to ~/TFLOW-0.9.1/tflow.sh test -t Test_Pipe")

    3. Then, for each component that can't be found:
        A. Install that component (if not already installed)
        B. Tell TFLOW where it is by editing TFLOW_Local_Settings.dat: 
           (Ex: vi ~/TFLOW-0.9.1/TFLOW_Local_Settings.dat)
               set [UTILITY]_LOCATION variable with executable conatining folder
           OR: set [UTILITY]_EXEC variable of full executable path
        C. When every component (that you care about) reports "Found!" you're ready for a project!

  - Start a Project
    1. Create a project directory
       Ex: mkdir ~/My_Project

    2. Copy desired pipe options file template to project directory:
       Ex: cp ~/TFLOW-0.9.1/Templates/options.dat.Trinity_Pipe_Template ~/My_Project/options.dat

    3. Set mandatory and suggested settings in options file:
       Ex: vi ~/My_Project/options.dat
           is_paired_reads <True/False>
           label           <tissue_name>
           BUSCO_type      <vertebrata>

    4. Run TFLOW:
        Ex: ~/TFLOW-0.9.1/tflow.sh run --raw_reads ~/My_Data/*.gz"
        "[TFLOW_DIR]/tflow.sh run --raw_left_reads ~/My_Data/*left_reads*.gz --raw_right_reads ~/My_Data/*right_reads*.gz"
    	    (template job submission scripts are included for each pipe type)
	 
    5. Track the Running Project:
        If you started TFLOW in the background (with "&" at the end),
            see what's going on in the job by running TFLOW in track mode.
            Ex:  cd ~/My_Project/
                 ~/TFLOW-0.9.1/tflow.sh track

  - Help

    For help, use the "-h" or "--help" flags.
        Ex: ~/TFLOW-0.9.1/tflow.sh run -h
        OR: ~/TFLOW-0.9.1/tflow.sh run --help


  - Advanced Use 
  
    For convience, you can:
        Add TFLOW dir to your system path with: "export PATH=[TFLOW_DIR]:$PATH"  
        Add TFLOW dir to your python path with: "export PYTHONPATH=[TFLOW_DIR]:$PYTHONPATH"
            (Replacing [TFLOW_DIR] with absolute path to TFLOW Directory Location)


 --- Operation ---
 
TFLOW 0.9.2 has 9 modes:

    run	    Run Job/Pipe
    track	    Track Job/Pipe
    analyze	    Analyze Job/Pipe
    read 	    Read Output of Job/Pipe
    test	    Test Job/Pipe
    stop	    Stop Running Job/Pipe
    clean	    Cleanup (Most) Standard Job/Pipe Files (Beta)
    reset	    Reset (Almost) All Job/Pipe Files Including Output (Beta)
    settings    Print Current TFLOW Settings and Exit

TFLOW 0.9.2 has 5 Supported Pipes:

    Trinity_Pipe         Trim Reads, Trinity-Assemble, CAP3-Assemble, CEGMA Analysis, 
                         and BUSCO Analysis
    Trimmed_Trinity_Pipe Trinity-Assemble, CAP3-Assemble, CEGMA Analysis, and BUSCO Analysis
    CAP3_Pipe	     CAP3-Assemble, CEGMA Analysis, and BUSCO Analysis
    Analysis_Pipe        Stat_Analysis, CEGMA Analysis and BUSCO Analysis
    Test_Pipe	     Non-Functional Pipe for Testing All Supported Segments

TFLOW 0.9.2 has 8 Supported Pipes Segments:

    Make_Read_Lists (v0.9)   Simple Parser From Read Files to Read File Lists
    Trimmomatic (v0.32)      Read Trimming Utilityad Output of Job/Pipe
    Trinity (v20140717)      De-Novo Transcriptome Assembler
    CAP3                     Sequence Assembler
    Package (v0.9)  	 Copy and Zip Final Sequence File Result
    Stat_Analysis (v0.9)     Perform Statistical Analysis on FASTA File Sequences
    CEGMA_Analysis (v2.2.29) Core Eukaryotic Gene Recapture Analysis via NCBI Blast
    BUSCO_Analysis (v2.2.29) Benchmarking Gene Recapture Analysis via NCBI Blast
    Summary (v0.9)  	 Summarize Pipe Results

 --- External Components ---
 
Components used for Pipe Segments. Newer versions of components are ideally also supported.

    - Trimmomatic -
    Used for trimming of raw sequence reads based on sequence quality and similarity
        to adapter sequences.
    Version:  0.3.2
    URL:      http://www.usadellab.org/cms/?page=trimmomatic'
    Citation: Bolger, Anthony M., Marc Lohse, and Bjoern Usadel. "Trimmomatic: a flexible 
              trimmer for Illumina sequence data." Bioinformatics (2014): btu170.

    - Trinity -
    De Novo Transcriptome Assembler, Assembles Raw Sequence Reads into Transcript Sequences.
    Version:  trinityrnaseq_r20140717
    URL:      http://trinityrnaseq.github.io/
    Citation: Full-length transcriptome assembly from RNA-Seq data without a reference genome.
              Grabherr MG, Haas BJ, Yassour M, Levin JZ, Thompson DA, Amit I, Adiconis X, Fan L,
	          Raychowdhury R, Zeng Q, Chen Z, Mauceli E, Hacohen N, Gnirke A, Rhind N, di Palma F,
	          Birren BW, Nusbaum C, Lindblad-Toh K, Friedman N, Regev A.
	          Nature Biotechnology 29, 644–652 (2011)
	          Paper: http://www.nature.com/nbt/journal/v29/n7/full/nbt.1883.html
	          Code:  http://trinityrnaseq.sf.net

    - CAP3 -
    DNA Sequence Assembler
    Version:  N/A
    URL:      http://seq.cs.iastate.edu/cap3.html
    Citation: Huang, X. and Madan, A. (1999) CAP3: A DNA sequence assembly program. 
    	  Genome Res., 9, 868-877.

    - NCBI BLAST -
    Sequence Comparison Tool
    Version:  2.2.29+
    URL:      http://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download
    Citation: Camacho C., Coulouris G., Avagyan V., Ma N., Papadopoulos J., Bealer K., 
              & Madden T.L. (2008) "BLAST+: architecture and applications." 
              BMC Bioinformatics 10:421


 --- Benchmarking Datasets ---
 
Datasets used for Gene Recapture Analysis.

    - CEGMA Sequence Database -
    Core Eukaryotic Genes Mapping Approach (CEGMA): Core Eukaryotic Genes Dataset
    Version:  Acquired 2015-04-22
    URL:      http://korflab.ucdavis.edu/datasets/cegma/
    Citation: Genis Parra, Keith Bradnam and Ian Korf. CEGMA: a pipeline to accurately annotate 
              core genes in eukaryotic genomes." Bioinformatics, 23: 1061-1067 (2007)

    - BUSCO Sequence Databases -
    Benchmarking sets of Universal Single-Copy Orthologs (BUSCO): Metazoa, Arthropoda, 
        Vertebrata, and Fungi Datasets
    Version:  OrthoDB7, Acquired 2015-04-22  
    URL:      ftp://cegg.unige.ch/OrthoDB7/BUSCO/
    Citation: Waterhouse et al, Nucleic Acids Research, January 2013, PMID:23180791
              OrthoDB: a hierarchical catalog of animal, fungal and bacterial orthologs.


 --- Options ---
 
The priority hierarchy of options is given by:
    Pipe Settings File > Command Line > Options File > Global Defaults > Segment Defaults  


--- Version History ---

2015-??-??, V.0.9.2

    Increased README Readability and added setup command examples.
    Added Termination of Tracking on Process Failure
    Added Internal Zipping and Unzipping of ".gz" Compressed Files
    Added "details" Mode to fasta_manip.py for Shortest/Longest Read Identification
    Added CSV outputting to Summary module
    Added prototype "reset" Mode for Deleting All Identified Segment Output Files
    Added prototype "clean" Mode for Deleting Unnecessary Segment Output Files
    Added "stop" Mode for Stopping Background Running Pipes, Segments, and External Components
    Bug Fixes

2015-05-12, V.0.9.1

    Added Stat_Analysis Segment
    Added Package Segment
    Added Test_Pipe Testing Pipeline
    Separated Pipes and Segments Into tflow/pipes and tflow/segments, Respectively
    Added Unpaired Read Trimming Mode
    Added Significantly More Code Commenting to Segments: BUSCO_Analysis, and CEGMA_Analysis
    Added Summary Module
    Added Warning When User-Supplied Setting is Overridden by Pipe Setting
    Fixed Trinity_Pipe Settings that Overrode BUSCO Tissue Type Selection
    Added test_TFLOW_setup.sh Convenience Script
    Added Significantly More Detail to README.txt
    Added CAP3_Pipe Job Templates
    Added CAP3 Assembly Pipe: CAP3_Pipe
    Minor Bug Fixes

2015-04-20, V.0.9.0

    Initial Beta Release
