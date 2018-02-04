#!/usr/bin/env python
import re, argparse, os
from recyclelib.utils import readfq

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'make_fasta_from_fastg converts fastg assembly graph to fasta format'
        )
    parser.add_argument('-g','--graph',
     help='(spades 3.50+) FASTG file to process [recommended: before_rr.fastg]',
     required=True, type=str
     )
    parser.add_argument('-o','--output',
     help='output file name for FASTA of cycles',
     required=False, type=str
     )

    return parser.parse_args()

def parse_lines(fastg, ofile):
    fp = open(fastg, 'r')
    for name,seq,qual in readfq(fp):
        name = re.sub('[:,]'," ", name[:-1]).split(" ")[0]
        if name[-1] == "'": continue # only need one of forward and reverse complement
        line = ">"+name+"\n"+seq+"\n"
        ofile.write(line)

if __name__=='__main__':
    args = parse_user_input()
    fastg = args.graph
    fp = open(fastg, 'r')
    files_dir = os.path.dirname(fp.name)

    # output 1 - fasta of sequences
    if args.output:
        fasta_ofile = args.output
    else:
        (root,ext) = os.path.splitext(fp.name)
        fasta_ofile = root + ext.replace(".fastg", ".nodes.fasta")

    f_nodes_fasta = open(fasta_ofile, 'w')
    parse_lines(fastg, f_nodes_fasta)
