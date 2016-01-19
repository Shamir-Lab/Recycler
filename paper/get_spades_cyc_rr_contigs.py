import argparse, sys, re, os
sys.path.insert(0, '../recycle/')
from recycle.utils import *
import numpy as np


def parse_user_input():
    parser = argparse.ArgumentParser(description='gets cycles from spades contigs output. Outputs these to separate fasta file')
    parser.add_argument('-g','--graph', help='Input (SPAdes 3.60+) FASTG to process', required=True, type=str)
    parser.add_argument('-p','--paths', help='Input (SPAdes 3.60+) contigs.paths file to process', required=True, type=str)
    parser.add_argument('-s','--sequences', help='Input (SPAdes 3.60+) contigs.fasta sequences to process', required=True, type=str)
    parser.add_argument('-m','--length', help='Minimum cycle length to keep (shorter cycles put in new graph file; default = 1000)', 
        required=False, type=int, default=1000)
    parser.add_argument('-k','--max_k', help='integer reflecting maximum k value used by the assembler',
     required=True, type=int, default=55)
    return parser.parse_args()

def get_cyc_labels_from_paths_file(paths):
	""" scans lines of paths file
		when length greater than 2 and 
		first node same as last, add label of path
		to returned set; add only forward strand versions
	"""
	labels = set({})
	fp = open(paths, 'r')
	lines = fp.read()
	fp.close()
	lines = re.sub(';\n',",", lines)
	lines = lines.split("\n")
	last_label = ''
	for ind in range(len(lines)):
		line = lines[ind].rstrip()
		if (ind % 4 == 2 or ind % 4 == 3):
			continue
		if ind%4==0:
			last_label = line
		else:
			path = line.split(",")
			if len(path)>1 and path[0]==path[-1]:
				labels.add(last_label)
	return labels	

def update_paths_dict(paths_dict, labels, seqs):
	fs = open(seqs, 'r')
	for name,seq,qual in readfq(fs):
		if name in labels:
			paths_dict[name] = seq
	return seqs



args = parse_user_input()
fastg = args.graph
paths = args.paths
seqs = args.sequences
min_length = args.length
max_k = args.max_k
fg = open(fastg, 'r')
files_dir = os.path.dirname(fg.name)

G = get_fastg_digraph(fastg)
SEQS = get_fastg_seqs_dict(fastg,G)
long_self_loops = get_long_self_loops(G, min_length, SEQS)
final_paths_dict = {}
path_count = 0

for nd in long_self_loops:
    name = get_spades_type_name(path_count, nd,
        SEQS, G, get_cov_from_spades_name(nd[0]))
    final_paths_dict[name] = nd
    path_count += 1

path_labels_set = get_cyc_labels_from_paths_file(paths)
update_paths_dict(final_paths_dict, path_labels_set, seqs)

# output - fasta of cyc sequences
fp = open(paths, 'r')
(root,ext) = os.path.splitext(fp.name)
fasta_ofile = root + ext.replace(".paths", ".cycs.fasta")
f_cycs_fasta = open(fasta_ofile, 'w')

for p in final_paths_dict.keys():
    if p not in path_labels_set:
    	# cycle = False b/c want last k-mer to match first
    	seq = get_seq_from_path(final_paths_dict[p], SEQS, max_k_val=max_k, cycle=False)
    else: 
    	seq = final_paths_dict[p]
    # 
    # print ""
    if len(seq)>=min_length:
        f_cycs_fasta.write(">" + p + "\n" + seq + "\n")
