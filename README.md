# Recycler

Recycler: an algorithm for detecting plasmids from de novo assembly graphs

##usage: recycle.py [-h] -i INPUT.FASTG -c COMP.FASTA [-l LENGTH] [-m MAX_CV]


##required arguments:
  -i, --input.fastg
                      Input (SPAdes 3.50+) FASTG to process (before_rr.fastg recommended)
  
  -c, --comp.fasta
                      FASTA of graph sequences to process [can be whole graph -- i.e., before_rr.fasta -- a single component, etc.]
                      
##optional arguments:
  -h, --help            show help message and exit
  
  -l, --length 
                        minimum length required for reporting [default: 1000]
  
  -m, --max_CV
                        coefficient of variation used for pre-selection
                        [default: 0.50, higher--> less restrictive]; Note: not
                        a requisite for selection
