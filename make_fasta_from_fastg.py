import re, argparse, os
from recycle.utils import readfq

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'recycle extracts cycles likely to be plasmids from metagenome and genome assembly graphs'
        )
    parser.add_argument('-g','--graph',
     help='(spades 3.50+) FASTG file to process [recommended: before_rr.fastg]',
     required=True, type=str
     )
    return parser.parse_args()

def parse_lines(fastg, ofile):
    lines = []
    fp = open(fastg, 'r')
    count = 0
    for name,seq,qual in readfq(fp):
        count += 1
        if count % 2 == 0: continue 
        name = re.sub('[:,]'," ", name[:-1]).split(" ")[0]
        line = ">"+name+"\n"+seq+"\n"
        ofile.write(line)
    
args = parse_user_input()
fastg = args.graph
fp = open(fastg, 'r')
files_dir = os.path.dirname(fp.name)

# output 1 - fasta of sequences
(root,ext) = os.path.splitext(fp.name)
fasta_ofile = root + ext.replace(".fastg", ".nodes.fasta")
f_nodes_fasta = open(fasta_ofile, 'w')
parse_lines(fastg, f_nodes_fasta)
