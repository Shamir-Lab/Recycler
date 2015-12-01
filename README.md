# Recycler

Recycler: an algorithm for detecting plasmids from de novo assembly graphs

##usage: python recycle.py [-h] -g GRAPH -s SEQUENCES [-l LENGTH] [-m MAX_CV]


##required arguments:
  -g, --GRAPH
                      Input FASTG to process (SPAdes 3.50 before_rr.fastg recommended)
  
  -s, --SEQUENCES
                      FASTA of graph sequences to process (can be whole graph or selected component)
                      
##optional arguments:
  -h, --help            show help message and exit
  
  -l, --length 
                        minimum length required for reporting [default: 1000]
  
  -m, --max_CV
                        coefficient of variation used for pre-selection
                        [default: 0.25, higher--> less restrictive]
