# Recycler

Recycler: an algorithm for detecting plasmids from de novo assembly graphs

##usage: python recycle.py [-h] -g GRAPH -s SEQUENCES [-l LENGTH] [-m MAX_CV]


##required arguments:
  -g, --GRAPH
                      Input FASTG to process (SPAdes 3.50 before_rr.fastg recommended)
  
  -s, --SEQUENCES
                      FASTA of graph sequences to process (can be whole graph or selected component); 
                      *Note:* node names must match those of graph file: e.g., **if you use before_rr.fastg as the graph,
                      please use before_rr.fasta as the sequence file** or a selection of sequences from it.
                      
##optional arguments:
  -h, --help            show help message and exit
  
  -l, --length 
                        minimum length required for reporting [default: 1000]
  
  -m, --max_CV
                        coefficient of variation used for pre-selection
                        [default: 0.25, higher--> less restrictive]

###outputs: 
\<prefix\>.cycs.fasta  - a fasta formatted file of predicted plasmids

\<prefix\>.cycs.paths_w_cov.txt - a text file containing information about each predicted plasmid.

The format for the second file is:
* *\<plasmid name\>* - e.g., RNODE_5_length_42666_cov_19.93685
* *\<node names in the original graph making up this cycle\>* - e.g., \('NODE_2801_length_42596_cov_19.8677', "NODE_2387_length_125_cov_34.7286'"\).
* *\<coverage levels of nodes at the time they are removed\>* - e.g., \[19.8677, 34.7286\]
* *\<node numbers in the original graph making up this cycle\>* - e.g., \[2801, 2387\]. This can be useful for visualizing the path in tools like [Bandage](https://rrwick.github.io/Bandage/)

