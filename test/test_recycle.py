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
	#### entry tested:
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
	# test removing again doesn't do anything
	update_node_coverage(G, test_node, 0)
	assert_equal(len(G.nodes()), num_nodes-2)
	assert_equal(len(G.edges()), num_edges-6)

def test_path_functions():
	# load test graph
	fastg = ROOT_DIR + "test/assembly_graph.fastg"
	test_node = "EDGE_1243_length_1496_cov_78.6919"
	G = get_fastg_digraph(fastg)
	comps = nx.strongly_connected_component_subgraphs(G)
	COMP = nx.DiGraph()

	# choose desired SCC based on node in it
	for c in comps:
		if test_node in c.nodes():
			COMP = c.copy()
			break
	SEQS = get_fastg_seqs_dict(fastg, COMP)
	# print COMP.nodes()
	# check sequences and nodes have been fetched
	assert_equal(len(COMP.nodes()), 8) # 3 cycle comp with isolated nodes removed
	assert_equal(len(SEQS), 8)
	assert_equal(SEQS["EDGE_675_length_69_cov_24.9286'"], 
		"TGTCCCTTTTACTGTTACAAAATGTCCCTTTTACTGTTACAAAATGTCCCTTTTACTGTTACAAAATGT")
	
	test_path = ('EDGE_1148_length_2822_cov_34.1811',
		'EDGE_71_length_961_cov_29.7759',
		'EDGE_1243_length_1496_cov_78.6919',
		"EDGE_69_length_2131_cov_28.8675'"
		)
	#69-, 71+, 1148+, 1243+, 676-, 677+, 1244- 
	test_fig8_path = ("EDGE_69_length_2131_cov_28.8675'",
		'EDGE_71_length_961_cov_29.7759',  'EDGE_1148_length_2822_cov_34.1811',
		'EDGE_1243_length_1496_cov_78.6919',
		"EDGE_676_length_1278_cov_32.0638'", 'EDGE_677_length_63_cov_57.625',
		 "EDGE_1244_length_5010_cov_35.8545'"
		)

	# test linear path, cycle, figure 8 all have correct length
	assert_equal(len(get_seq_from_path(test_path, SEQS, cycle=True)), 7190)	
	assert_equal(len(get_seq_from_path(test_path, SEQS, cycle=False)), 7245)
	assert_equal(len(get_seq_from_path(test_fig8_path, SEQS, cycle=True)), 13376)
	assert_equal(len(get_seq_from_path(test_fig8_path, SEQS, cycle=False)), 13431)

	# single node path - CV = 0
	assert_equal(get_wgtd_path_coverage_CV(("EDGE_69_length_2131_cov_28.8675'",),
		G, SEQS),0)

	# set coverage values to dummy values for testing
	for node in test_path:
		update_node_coverage(G, node, 10)
	# const coverage path - CV = 0
	assert_equal(get_wgtd_path_coverage_CV(test_path, G, SEQS), 0)
	# CV doesn't depend on magnitude
	COMP2 = COMP.copy()
	COMP3 = COMP.copy()

	for node in test_path:
		update_node_coverage(G, node, get_cov_from_spades_name(node))
		update_node_coverage(COMP2, node, get_cov_from_spades_name(node)*10)
		update_node_coverage(COMP3, node, get_cov_from_spades_name(node)*1000)
	
	assert_equal(get_wgtd_path_coverage_CV(test_path, COMP2, SEQS), 
		get_wgtd_path_coverage_CV(test_path, COMP3, SEQS))

	assert_equal(get_total_path_mass(test_path, COMP3) 
		/ get_total_path_mass(test_path, COMP2), 100)


def test_get_long_self_loops():
	fastg = ROOT_DIR + "test/assembly_graph.fastg"
	G = get_fastg_digraph(fastg)
	SEQS = get_fastg_seqs_dict(fastg, G)
	min_length = 1000 
	# test returned set
	assert_true(("EDGE_2131_length_56011_cov_21.811'",) in get_long_self_loops(G,min_length,SEQS))
	assert_true(("EDGE_299_length_56_cov_728",) not in get_long_self_loops(G,min_length,SEQS)) # too short
	assert_true(("EDGE_1548_length_7806_cov_2.3197'",) not in get_long_self_loops(G,min_length,SEQS)) # linear

	# test loop nodes removed from graph
	assert_true("EDGE_2131_length_56011_cov_21.811'" not in G)
	assert_true("EDGE_2131_length_56011_cov_21.811" not in G)
	assert_true("EDGE_299_length_56_cov_728" not in G)
	assert_true("EDGE_299_length_56_cov_728'" not in G)
