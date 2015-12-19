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

def test_get_fastg_digraph():
	fastg = ROOT_DIR + "test/assembly_graph.fastg"
	G = get_fastg_digraph(fastg)
	assert_equal(len(G.nodes()), 2*3331)
	assert_equal(len(G.edges()), 2*905)

# get_fastg_digraph(