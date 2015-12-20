from nose.tools import *
from recycle.utils import *

ROOT_DIR = "/specific/a/home/cc/cs/rozovr/recycle/"

def test_rc():
	assert_equal("ACGTT", rc_seq("AACGT"))
	assert_equal("ACGT", rc_seq("ACGT")) # palindrome
	assert_false("ACGT" == rc_seq("AAAA"))

def test_spades_name_functions():
	test_name = "RNODE_19_length_15528_cov_23.37027"
	num, length, cov = (get_num_from_spades_name(test_name), 
		get_length_from_spades_name(test_name), 
		get_cov_from_spades_name(test_name))
	assert_equal((num,length,cov), (19,15528,23.37027))
	assert_equal("RNODE_19_length_15528_cov_23.37027'", rc_node(test_name))
	assert_equal(test_name, rc_node("RNODE_19_length_15528_cov_23.37027'"))


def test_get_fastg_digraph():
	fastg = ROOT_DIR + "test/assembly_graph.fastg"
	G = get_fastg_digraph(fastg)
	assert_equal(len(G.nodes()), 2*3331)
	assert_equal(len(G.edges()), 2*905)

def test_coverage_funcs():
	# entry tested:
	# >EDGE_184_length_56_cov_59:EDGE_1102_length_76_cov_138.143';
	# TGGGTCGTTAGCTCAGTTGGTAGAGCAGTTGACTTTTAATCAATTGGTCGCAGGTT
	# >EDGE_184_length_56_cov_59':EDGE_182_length_155_cov_24.71,EDGE_183_length_3548_cov_23.29;
	# AACCTGCGACCAATTGATTAAAAGTCAACTGCTCTACCAACTGAGCTAACGACCCA

	# load sample graph, store initial sizes
	fastg = ROOT_DIR + "test/assembly_graph.fastg"
	G = get_fastg_digraph(fastg)
	test_node = 'EDGE_184_length_56_cov_59'
	test_node_rc = rc_node(test_node)
	num_nodes = len(G.nodes())
	num_edges = len(G.edges())

	# initial coverage is spades coverage
	assert_equal(get_cov_from_spades_name_and_graph(test_node, G), 59)

	# update coverage value, test again
	update_node_coverage(G, test_node, 10)
	assert_equal(get_cov_from_spades_name_and_graph(test_node, G), 10)
	assert_equal(get_cov_from_spades_name_and_graph(test_node_rc, G), 10)

	# remove node by updating coverage to 0
	# test both F and R are removed, test all edges to 
	# node are removed
	update_node_coverage(G, test_node, 0)
	assert_equal(get_cov_from_spades_name_and_graph(test_node, G), 0)
	assert_equal(get_cov_from_spades_name_and_graph(test_node_rc, G), 0)
	# F and R versions of both nodes and edges removed
	assert_equal(len(G.nodes()), num_nodes-2)
	assert_equal(len(G.edges()), num_edges-6)


# def get_cov_from_spades_name_and_graph(name,G):
#     if name not in G:
#         return 0
#     if 'cov' in G.node[name]:
#         return G.node[name]['cov']
#     else:
#         return get_cov_from_spades_name(name)

# def get_spades_base_mass(G, name):
#     length = get_length_from_spades_name(name)
#     coverage = get_cov_from_spades_name_and_graph(name,G)
#     return length * coverage

# def get_seq_from_path(path, seqs, max_k_val=55):
#     seq = get_fasta_stranded_seq(seqs, path[0])
#     if len(path)!=1:
#         for p in path[1:]:
#             seq += get_fasta_stranded_seq(seqs, p)[max_k_val:]
#     return seq

# def get_total_path_length(path, seqs):
#     # return sum([get_length_from_spades_name(n) for n in path])
#     seq = get_seq_from_path(path, seqs)
#     return len(seq)

# def get_wgtd_path_coverage_CV(path, G, seqs, max_k_val=55):
#     covs = np.array([get_cov_from_spades_name_and_graph(n,G) for n in path])
#     if len(covs)< 2: return 0.000001
#     # mean = np.mean(covs)
#     wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
#     tot_len = get_total_path_length(path, seqs)
#     wgts = np.multiply(wgts, 1./tot_len)
#     mean = np.average(covs, weights = wgts)
#     # try:
#     # diffs = covs - mean    
#     std = np.sqrt(np.dot(wgts,(covs-mean)**2))
#     # except TypeError:
#     #     print type(wgts), type(diffs**2), wgts, covs-mean
#     return std/mean


# def get_path_coverage_CV(path,G):
#     covs = np.array([get_cov_from_spades_name_and_graph(n,G) for n in path])
#     if len(covs)< 2: return 0.000001
#     mean = np.mean(covs)
#     std = np.std(covs)
#     # if mean == 0: return 1000 
#     return std/mean

# def get_total_path_mass(path,G):
#     return sum([get_length_from_spades_name(p) * \
#         get_cov_from_spades_name_and_graph(p,G) for p in path])


# # get_fastg_digraph(