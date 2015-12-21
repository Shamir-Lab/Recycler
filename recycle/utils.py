import numpy as np
import networkx as nx
import re
complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}


def readfq(fp): # this is a generator function
    """ # lh3's fast fastX reader: 
        https://github.com/lh3/readfq/blob/master/readfq.py
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def rc_seq(dna):
    rev = reversed(dna)
    return "".join([complements[i] for i in rev])

def get_num_from_spades_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[1]
    return int(contig_length)      

def get_length_from_spades_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[3]
    return int(contig_length)

def get_cov_from_spades_name(name):
    name_parts = name.split("_")
    cov = name_parts[5]
    if cov[-1]=="'": cov=cov[:-1]
    return float(cov)

def get_fastg_digraph(fastg_name):
    """ scans through fastg headers as an adjacency list
        builds and returns a nx directed graph using adjacencies
        note: no connections are created between each node and its 
        rc node - we need to take care to maintain these 
    """
    lines = []
    fp = open(fastg_name, 'r')
    for name,seq,qual in readfq(fp):
        name = re.sub('[:,]'," ", name[:-1])
        lines.append(name)
    G = nx.DiGraph()
    return nx.parse_adjlist(lines, create_using=G)

def get_fastg_seqs_dict(fastg_name, G):
    """ returns a dictionary of sequences in graph 
        where node names are keys and sequence strings
        are values; useful for saving memory when G
        is a subgraph (e.g., a component)
    """ 
    fp = open(fastg_name, 'r')
    seqs = {}
    for name,seq,qual in readfq(fp):
        name_parts = re.sub('[:,]'," ", name[:-1]).split()
        node = name_parts[0]
        if node in G.nodes():
            seqs[node] = seq
    return seqs

def rc_node(node):
    """ gets reverse complement
        spades node label
    """ 
    if node[-1] == "'": return node[:-1]
    else: return node + "'"


def get_cov_from_spades_name_and_graph(name,G):
    if name not in G:
        return 0
    if 'cov' in G.node[name]:
        return G.node[name]['cov']
    else:
        return get_cov_from_spades_name(name)

def update_node_coverage(G, node, new_cov):
    """ changes coverage value stored in 'cov'
        field on both F and R version of a node
        if new_cov is 0, both node versions are removed
    """
    if node not in G.nodes(): # nothing to be done, perhaps already removed
        return
    if new_cov == 0: 
        G.remove_node(node)
        G.remove_node(rc_node(node))
    else:
        G.add_node(node, cov=new_cov)
        G.add_node(rc_node(node), cov=new_cov)

def get_spades_base_mass(G, name):
    length = get_length_from_spades_name(name)
    coverage = get_cov_from_spades_name_and_graph(name,G)
    return length * coverage

def get_seq_from_path(path, seqs, max_k_val=55, cycle=True):
    """ retrieves sequence from a path;
        instead of specifying cycles by having the first and 
        last node be equal, the user must decide if a cycle is
        the intent to avoid redundant k-mers at the ends
    """
    start = seqs[path[0]]
    if len(path)==1:
        if cycle: 
            return start[max_k_val:]
        else:
            return start
    else:
        seq = ''
        for p in path:
            seq += seqs[p][max_k_val:]
        if cycle: return seq
        else: return start[:max_k_val] + seq

# def get_total_path_length(path, seqs):
#     # return sum([get_length_from_spades_name(n) for n in path])
#     seq = get_seq_from_path(path, seqs)
#     return len(seq)

def get_wgtd_path_coverage_CV(path, G, seqs, max_k_val=55):
    if len(path)< 2: return 0
    covs = np.array([get_cov_from_spades_name_and_graph(n,G) for n in path])
    wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
    tot_len = len(get_seq_from_path(path, seqs, cycle=True))
    wgts = np.multiply(wgts, 1./tot_len)
    mean = np.average(covs, weights = wgts)
      
    std = np.sqrt(np.dot(wgts,(covs-mean)**2))
    return std/mean


def get_path_coverage_CV(path,G):
    covs = np.array([get_cov_from_spades_name_and_graph(n,G) for n in path])
    if len(covs)< 2: return 0.000001
    mean = np.mean(covs)
    std = np.std(covs)
    # if mean == 0: return 1000 
    return std/mean

def get_total_path_mass(path,G):
    return sum([get_length_from_spades_name(p) * \
        get_cov_from_spades_name_and_graph(p,G) for p in path])

# def get_fasta_stranded_seq(seqs, seq_name):
#     """ gets sequence corresponding 
#         to same strand as fasta input file 
#         or rc seq depending on sequence name
#     """
#     if seq_name[-1]!="'":
#         return seqs[seq_name]
#     else: 
#         return rc_seq(seqs[seq_name[:-1]])