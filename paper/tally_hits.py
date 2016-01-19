import argparse, sys
sys.path.insert(0, '../recycle/')
from recycle.utils import *
import numpy as np


def count_property_range_hits(prop, node_dict, hits):
	""" picks which values to use in tuples based on property
		counts vals having min_val <= val < max_val
	 	unless max_val == overall_max, where min_val <= val <= max_val 
	 	is used instead
	"""
	res = []
	# sets tuple position to use in dict value
	switcher = {
        "length": (0,(0,4000,8000,12000,16000,20000)),
        "steps": (1,(0,2,4,8,16,32)),
        "cov": (2,(1,10,100,1000,10000,100000)),
        "cv": (3, (0,0.05,0.10,0.15,0.20,0.25))
    }
	if prop not in switcher:
		return res
	tup_pos = switcher[prop][0]
	node_cnt = 0
	pos_cnt = 0
	for ind in range(len(switcher[prop][1])-1):
		min_val = switcher[prop][1][ind]
		max_val = switcher[prop][1][ind+1]
		for node in node_dict.keys():
			val = node_dict[node][tup_pos]
			if ind < len(switcher[prop][1])-2:
				range_test_val = (min_val <= val < max_val)
			else:
				range_test_val = (min_val <= val <= max_val)
			# print "range bool is", range_test_val
			if range_test_val:
				node_cnt += 1
				if node in hits: pos_cnt += 1
		if node_cnt > 0:
			res.append( (pos_cnt, node_cnt, round(float(pos_cnt)/node_cnt,2)))
		else:
			res.append((0,0,0))
		node_cnt = 0
		pos_cnt = 0
	return res

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'counts number of candidates vs TP hits for a certain property and range combination'
        )
    parser.add_argument('-p','--pref',
     help='prefix to recycler outputs',
     required=True, type=str
     )
    parser.add_argument('-n',
        '--nodes', help='nodes list including accepted hits to reference sequences',
         required=True, type=str
         )
    
    return parser.parse_args()


def get_path_vals_cv(path, covs, max_k_val=55):
    """ returns cv value based on coverage values
    	when path was removed - stored in 
    	paths_w_cov.txt output file
    """
    wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
    tot_len = sum(wgts)
    wgts = np.multiply(wgts, 1./tot_len)
    mean = np.average(covs, weights = wgts)
    # try:
    # diffs = covs - mean    
    std = np.sqrt(np.dot(wgts,(covs-mean)**2))
    return std/mean


############### ENTRY POINT ####################
# inputs: paths_w_cov.txt file, 
args = parse_user_input()

paths_file = args.pref + '.cycs.paths_w_cov.txt'

# create dict of RNODES as keys, values
# as tuples of (total_length, num_steps, coverage)
f = open(paths_file, 'r')
lines = f.readlines()
rnode_dict = {}
for ind in range(len(lines)/4):
	name = lines[ind*4].rstrip()
	path = lines[ind*4 + 1].rstrip()
	path_covs = np.array([float(a) for a in lines[ind*4 + 2].rstrip()[1:-1].split(",")])

	cov = get_cov_from_spades_name(name)
	length = get_length_from_spades_name(name)
	num_steps = len(path.split(','))
	cv = get_path_vals_cv(path[1:-1].split(","), path_covs)
	# print name, num_steps, path
	rnode_dict[name] = (length, num_steps, cov, cv)

# issue - single node path RNODEs are not in path_w_cov files
# read in cycs.fasta file, add back single node paths 
cycs_file = args.pref + '.cycs.fasta'
f = open(cycs_file, 'r')
lines = f.readlines()
for ind in range(len(lines)/2):
	name = lines[ind*2][1:].rstrip()
	if name in rnode_dict: continue
	else:
		cov = get_cov_from_spades_name(name)
		length = get_length_from_spades_name(name)
		# print name, 1
		rnode_dict[name] = (length, 1, cov, 0)


# nucmer.delta file parsed to RNODES having 100/80 hits with 
#/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l /home/nasheran/rozovr/recycle_paper_data/ref_800/before_rr.nucmer.delta | awk '$10==100.00 && $15>=80.00' | cut -d'|' --complement -f 1-6 | cut -f2 > /home/nasheran/rozovr/recycle_paper_data/rnode_hits.txt
# need to read names in and get set
hits_file = args.nodes
f = open(hits_file,'r')
hits = set([])
lines = f.readlines()
for line in lines:
	hits.add(line.rstrip())


print "length: ", count_property_range_hits("length", rnode_dict, hits)
print "steps: ", count_property_range_hits("steps", rnode_dict, hits)
print "coverage: ", count_property_range_hits("cov", rnode_dict, hits)
print "CV: ", count_property_range_hits("cv", rnode_dict, hits)


