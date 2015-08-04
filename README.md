# recycle
usage: recycle.py [-h] -i INPUT.FASTG -c COMP.FASTA [-l LENGTH] [-m MAX_CV]
<br>
recycle extracts cycles likely to be plasmids from metagenome and genome
assembly graphs
<br>
required arguments:<br>
  -i INPUT.FASTG, --input.fastg INPUT.FASTG<br>
                        Input (SPAdes 3.50+) FASTG to process
  -c COMP.FASTA, --comp.fasta COMP.FASTA<br>
                        Input graph component FASTA to process [can be whole
                        graph -- not recommended]
                        <br>
optional arguments:<br>
  -h, --help            show help message and exit<br>
  -l LENGTH, --length LENGTH<br>
                        minimum length required for reporting [default: 1000]
  -m MAX_CV, --max_CV MAX_CV<br>
                        coefficient of variation used for pre-selection
                        [default: 0.50, higher--> less restrictive]; Note: not
                        a requisite for selection
