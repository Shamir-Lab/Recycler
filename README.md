# recycle
usage: recycle.py [-h] -i INPUT.FASTG -c COMP.FASTA [-l LENGTH] [-m MAX_CV]

recycle extracts cycles likely to be plasmids from metagenome and genome
assembly graphs

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT.FASTG, --input.fastg INPUT.FASTG
                        Input (SPAdes 3.50+) FASTG to process
  -c COMP.FASTA, --comp.fasta COMP.FASTA
                        Input graph component FASTA to process [can be whole
                        graph -- not recommended]
  -l LENGTH, --length LENGTH
                        minimum length required for reporting [default: 1000]
  -m MAX_CV, --max_CV MAX_CV
                        coefficient of variation used for pre-selection
                        [default: 0.50, higher--> less restrictive]; Note: not
                        a requisite for selection
