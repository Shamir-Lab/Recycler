from utils import *

# print np.mean([1,2])

def get_mod_seq(dna, start, end, limit):
	if start >= limit and end >= limit:
		return dna[(start%limit):(end%limit)]
	elif end >= limit:
		return dna[start:limit] + dna[:(end%limit)]
	else:
		return dna[start:end]

def get_cyclic_kmers(seq, k):
	res = []
	for i in range(len(seq)):
		res.append(get_mod_seq(seq,i,i+k,len(seq)))
	return res
	# yield get_mod_seq(seq,i,i+k,len(seq))


def test_cyclic_kmers():
	assert(get_cyclic_kmers("ACGTC", 3)==["ACG", "CGT", "GTC", "TCA", "CAC"])

def canonicalize_kmers(kmers):
	res = []
	for k in kmers:
		res.append(min(k,rc_seq(k)))
	return res

def test_canonicalize_kmers():
	assert(canonicalize_kmers(["ACG", "CGT", "GTC", "TCA", "CAC"])==["ACG", "ACG", "GAC", "TCA", "CAC"])

def get_canonical_cyclic_seq_list(seq,k):
	""" returns a sorted list representation of
		canonical kmers in a cyclic sequence - used 
		for comparing between sequences when starting points
		and strand may differ 
	"""
	return canonicalize_kmers(get_cyclic_kmers(seq,k)).sort()

def test_get_canon_cycle_seq_list():
	# test rc case
	assert(get_canonical_cyclic_seq_list("ACGTC", 3)==get_canonical_cyclic_seq_list("GACGT", 3))
	# test rotated case
	assert(get_canonical_cyclic_seq_list("ACGTC", 3)==get_canonical_cyclic_seq_list("CACGT", 3))
	# test rotated rc case
	assert(get_canonical_cyclic_seq_list("ACGTC", 3)==get_canonical_cyclic_seq_list("TGACG", 3))

# read in reference seqs, rotate them
# to canonicalize



# fp = open(seqs_name, 'r')
# seq_nodes = []
# seqs = {}
# for name,seq,qual in readfq(fp):
#     # avoid sources & sinks
#     # TODO: get rid of higher order sources/sinks - e.g., 
#     # sinks caused by removal of sinks, ...
#     if G.out_degree(name)!=0 and G.in_degree(name)!=0:
#         seq_nodes.extend([name, name+"'"])
#         seqs[name] = seq

test_cyclic_kmers()
test_canonicalize_kmers()
test_get_canon_cycle_seq_list()
