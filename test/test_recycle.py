from nose.tools import *
from recyclelib.utils import *

ROOT_DIR = "test/"
TEST_PATH = ('EDGE_1148_length_2822_cov_34.1811',
		'EDGE_71_length_961_cov_29.7759',
		'EDGE_1243_length_1496_cov_78.6919',
		"EDGE_69_length_2131_cov_28.8675'"
		)
TEST_FIG8_PATH = ("EDGE_69_length_2131_cov_28.8675'",
		'EDGE_71_length_961_cov_29.7759',  'EDGE_1148_length_2822_cov_34.1811',
		'EDGE_1243_length_1496_cov_78.6919',
		"EDGE_676_length_1278_cov_32.0638'", 'EDGE_677_length_63_cov_57.625',
		 "EDGE_1244_length_5010_cov_35.8545'"
		)
# plasmid in E2022 requiring both + and - node traversal
# 286+, 287+, 468-, 469+, 476+, 286-
TEST_REPEAT_PATH = ('EDGE_286_length_92_cov_109.162',
 "EDGE_468_length_4093_cov_54.6159'", 'EDGE_469_length_420_cov_90.6849', 
 'EDGE_476_length_2009_cov_45.5921', "EDGE_286_length_92_cov_109.162'", 
 'EDGE_287_length_26768_cov_58.5435'
 )


def get_sample_graph_comp_seqs():
	# load test graph
	fastg = ROOT_DIR + "assembly_graph.fastg"
	test_node = "EDGE_1243_length_1496_cov_78.6919"
	G = get_fastg_digraph(fastg)
	comps = nx.strongly_connected_component_subgraphs(G)
	COMP = nx.DiGraph()

	# choose desired SCC based on node in it
	for c in comps:
		if test_node in c.nodes():
			COMP = c.copy()
			break
	SEQS = get_fastg_seqs_dict(fastg, G)
	return G,COMP,SEQS

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
	fastg = ROOT_DIR + "assembly_graph.fastg"
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
	fastg = ROOT_DIR + "assembly_graph.fastg"
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

def TEST_PATH_functions():
	G,COMP,SEQS = get_sample_graph_comp_seqs()
	# check sequences and nodes have been fetched
	assert_equal(len(COMP.nodes()), 8) # 3 cycle comp with isolated nodes removed
	assert_equal(SEQS["EDGE_675_length_69_cov_24.9286'"], 
		"TGTCCCTTTTACTGTTACAAAATGTCCCTTTTACTGTTACAAAATGTCCCTTTTACTGTTACAAAATGT")
	
	#69-, 71+, 1148+, 1243+, 676-, 677+, 1244- 
	

	# test linear path, cycle, figure 8 all have correct length
	assert_equal(len(get_seq_from_path(TEST_PATH, SEQS, cycle=True)), 7190)	
	assert_equal(len(get_seq_from_path(TEST_PATH, SEQS, cycle=False)), 7245)
	assert_equal(len(get_seq_from_path(TEST_FIG8_PATH, SEQS, cycle=True)), 13376)
	assert_equal(len(get_seq_from_path(TEST_FIG8_PATH, SEQS, cycle=False)), 13431)

	# single node path - CV = 0
	assert_equal(get_wgtd_path_coverage_CV(("EDGE_69_length_2131_cov_28.8675'",),
		G, SEQS),0)

	# set coverage values to dummy values for testing
	for node in TEST_PATH:
		update_node_coverage(G, node, 10)
	# const coverage path - CV = 0
	assert_equal(get_wgtd_path_coverage_CV(TEST_PATH, G, SEQS), 0)
	# CV doesn't depend on magnitude
	COMP2 = COMP.copy()
	COMP3 = COMP.copy()

	for node in TEST_PATH:
		update_node_coverage(G, node, get_cov_from_spades_name(node))
		update_node_coverage(COMP2, node, get_cov_from_spades_name(node)*10)
		update_node_coverage(COMP3, node, get_cov_from_spades_name(node)*1000)
	
	assert_equal(get_wgtd_path_coverage_CV(TEST_PATH, COMP2, SEQS), 
		get_wgtd_path_coverage_CV(TEST_PATH, COMP3, SEQS))

	assert_equal(get_total_path_mass(TEST_PATH, COMP3) 
		/ get_total_path_mass(TEST_PATH, COMP2), 100)


def test_get_long_self_loops():
	G,COMP,SEQS = get_sample_graph_comp_seqs()

	min_length = 1000 
	# test returned set
	assert_true(("EDGE_2131_length_56011_cov_21.811",) in get_long_self_loops(G,min_length,SEQS))
	assert_true(("EDGE_2131_length_56011_cov_21.811'",) not in get_long_self_loops(G,min_length,SEQS)) # rc version not in
	assert_true(("EDGE_299_length_56_cov_728",) not in get_long_self_loops(G,min_length,SEQS)) # too short
	assert_true(("EDGE_1548_length_7806_cov_2.3197'",) not in get_long_self_loops(G,min_length,SEQS)) # linear

	# test loop nodes removed from graph
	assert_true("EDGE_2131_length_56011_cov_21.811'" not in G)
	assert_true("EDGE_2131_length_56011_cov_21.811" not in G)
	assert_true("EDGE_299_length_56_cov_728" not in G)
	assert_true("EDGE_299_length_56_cov_728'" not in G)

def test_get_unoriented_sorted_str():
	TEST_PATH = ('EDGE_1148_length_2822_cov_34.1811',
		'EDGE_71_length_961_cov_29.7759',
		'EDGE_1243_length_1496_cov_78.6919',
		"EDGE_69_length_2131_cov_28.8675'"
		)
	assert_equal(get_unoriented_sorted_str(TEST_PATH), 
		"EDGE_1148_length_2822_cov_34.1811'EDGE_1243_length_1496_cov_78.6919'EDGE_69_length_2131_cov_28.8675'EDGE_71_length_961_cov_29.7759'")

def test_enum_high_mass_shorTEST_PATHs():
	# get component, gen. cycles on it
	## should refactor this out to function, 
	# as long as I know nosetests won't call that function
	G,COMP,SEQS = get_sample_graph_comp_seqs()

	paths = enum_high_mass_shortest_paths(COMP)
	# print paths
	assert_true(len(paths) < len(COMP.nodes()))

	
	for n in TEST_PATH:
		update_node_coverage(COMP, n, 0)

	paths = enum_high_mass_shortest_paths(COMP)
	assert_equal(len(paths), 1) # only 2 contig loop remains
	assert_true(paths[0] in (('EDGE_677_length_63_cov_57.625', "EDGE_675_length_69_cov_24.9286'"), 
		("EDGE_675_length_69_cov_24.9286'", 'EDGE_677_length_63_cov_57.625')))

def test_get_non_repeat_nodes():
	G,COMP,SEQS = get_sample_graph_comp_seqs()
	
	assert_true(get_non_repeat_nodes(G,TEST_PATH), ['EDGE_71_length_961_cov_29.7759'])

	sing_nodes = get_non_repeat_nodes(G,TEST_FIG8_PATH)
	assert_true('EDGE_71_length_961_cov_29.7759' in sing_nodes)
	assert_true("EDGE_676_length_1278_cov_32.0638'" in sing_nodes)
	assert_true("EDGE_1244_length_5010_cov_35.8545'" in sing_nodes)

def test_get_contigs_of_mates():
	G,COMP,SEQS = get_sample_graph_comp_seqs()
	bamfile = pysam.AlignmentFile(ROOT_DIR+"test.sort.bam", 'rb')
	# note mapped to positive nodes; thus only their names valid
	# also, need to change labels starting with "EDGE_" to "NODE_"
	mate_tigs = get_contigs_of_mates("EDGE_1244_length_5010_cov_35.8545", bamfile, G)
	assert_true("EDGE_676_length_1278_cov_32.0638" in mate_tigs)


def test_is_good_cyc():
	G,COMP,SEQS = get_sample_graph_comp_seqs()
	bamfile = pysam.AlignmentFile(ROOT_DIR+"test.sort.bam", 'rb')
	mate_tigs = get_contigs_of_mates("EDGE_800_length_15304_cov_22.6688", bamfile, G)

	assert_false(is_good_cyc(('EDGE_800_length_15304_cov_22.6688', "EDGE_801_length_279_cov_57.2411'"), G, bamfile))		
	assert_true(is_good_cyc(TEST_PATH, G, bamfile))	

def test_get_path_covs():
	fastg = ROOT_DIR+"assembly_graph2.fastg"
	G = get_fastg_digraph(fastg)
	assert_equal(get_node_cnts_hist(TEST_REPEAT_PATH)['EDGE_286_length_92_cov_109.162'], 2)
	assert_true(all(get_path_covs(TEST_REPEAT_PATH,G) == np.array([109.162/2, 54.6159, 90.6849, 45.5921, 109.162/2, 58.5435])))
	# ('EDGE_286_length_92_cov_109.162', "EDGE_468_length_4093_cov_54.6159'", 'EDGE_469_length_420_cov_90.6849', 'EDGE_476_length_2009_cov_45.5921', "EDGE_286_length_92_cov_109.162'", 'EDGE_287_length_26768_cov_58.5435')
