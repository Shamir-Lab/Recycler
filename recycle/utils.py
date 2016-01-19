import numpy as np
import networkx as nx
import re, pysam
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
        if rc_node(node) in G.nodes():
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

def get_wgtd_path_coverage_CV(path, G, seqs, max_k_val=55):
    if len(path)< 2: return 0
    mean, std = get_path_mean_std(path, G, seqs, max_k_val)
    return std/mean

def get_node_cnts_hist(path):
    d = {}
    for p in path:
        # always count based on positive node
        pos_name = p if (p[-1]!="'") else p[:-1]    
        d[pos_name] = d.get(pos_name,0) + 1
    return d

def get_path_covs(path,G):
    covs = [get_cov_from_spades_name_and_graph(n,G) for n in path]
    cnts = get_node_cnts_hist(path)
    for i in range(len(path)):
        p = path[i]
        pos_name = p if (p[-1]!="'") else p[:-1]
        if cnts[pos_name] > 1:
            covs[i] /= cnts[pos_name]
    return covs


def get_path_mean_std(path, G, seqs, max_k_val=55):
    # covs = np.array([get_cov_from_spades_name_and_graph(n,G) for n in path])
    covs = get_path_covs(path,G)
    wgts = np.array([(get_length_from_spades_name(n)-max_k_val) for n in path])
    tot_len = len(get_seq_from_path(path, seqs, cycle=True))
    wgts = np.multiply(wgts, 1./tot_len)
    mean = np.average(covs, weights = wgts)
    std = np.sqrt(np.dot(wgts,(covs-mean)**2))
    return (mean,std)

def update_path_coverage_vals(path, G, seqs):
    mean, _ = get_path_mean_std(path, G, seqs)
    # covs = np.array([get_cov_from_spades_name_and_graph(n,G) for n in path])
    covs = get_path_covs(path,G)
    new_covs = covs - mean
    for i in range(len(path)):
        if new_covs[i] > 0:
            update_node_coverage(G,path[i],new_covs[i])
        else:
            update_node_coverage(G,path[i],0)

def get_total_path_mass(path,G):
    return sum([get_length_from_spades_name(p) * \
        get_cov_from_spades_name_and_graph(p,G) for p in path])

def get_long_self_loops(G, min_length, seqs):
    """ returns set of self loop nodes paths that are longer 
        than min length; removes those and short self loops
        from G
    """
    self_loops = set([])
    to_remove = []

    for nd in G.nodes_with_selfloops():
        nd_path = (nd,)
        if len(get_seq_from_path(nd_path, seqs)) >= min_length \
        and (rc_node(nd),) not in self_loops:
            self_loops.add(nd_path)                
        to_remove.append(nd)

    for nd in to_remove:
        update_node_coverage(G, nd, 0)
    return self_loops

def get_unoriented_sorted_str(path):
    """ creates unique, orientation-oblivious string representation of path, 
        used to make sure node covered whenever rc of node is; 
        lets us avoid issue of rc of node having different weight than node
    """
    all_rc_path = []
    for p in path:
        if p[-1] != "'": p = p+"'"
        all_rc_path.append(p)
    return "".join(sorted(all_rc_path))

def enum_high_mass_shortest_paths(G, seen_paths=None):
    """ given component subgraph, returns list of paths that
        - is non-redundant (includes) no repeats of same cycle
        - includes all shortest paths starting at each node n (assigning 
        node weights to be 1/(length * coverage)) to each of 
        its predecessors, and returning to n
    """
    if seen_paths == None:
        seen_paths = []
    nodes = []
    nodes = list(G.nodes()) # creates a copy

    unq_sorted_paths = set([])
    # in case orientation obliv. sorted path strings passed in
    for p in seen_paths:
        unq_sorted_paths.add(p)
    paths = []
    # use add_edge to assign edge weights to be 1/mass of starting node
    for e in G.edges():
        G.add_edge(e[0], e[1], cost = 1./get_spades_base_mass(G, e[0]))

    for node in nodes:
        # if node[-1] == "'": continue
        for pred in G.predecessors(node):
            # needed because some nodes removed on updates
            if pred not in G: continue

            try:
                path = nx.shortest_path(G, source=node,
                    target=pred, weight='cost')
            except nx.exception.NetworkXNoPath:
                continue

            
            # below: create copy of path with each node as rc version
            # use as unique representation of a path and rc of its whole            
            unoriented_sorted_path_str = get_unoriented_sorted_str(path)

            # here we avoid considering cyclic rotations of identical paths
            # by sorting their string representations (all_rc_path above) 
            # and comparing against the set already stored
            if unoriented_sorted_path_str not in unq_sorted_paths:
                unq_sorted_paths.add(unoriented_sorted_path_str)
                paths.append(tuple(path))
    
    return paths

def get_non_repeat_nodes(G, path):
    """ returns a list of all non-repeat (in degree and out-degree
        == 1) nodes in a path; if there are no such nodes,
        returns an empty list
        NB: G input should be whole graph, not specific SCC, to avoid 
        disregarding isolated nodes  
    """
    sing_nodes = []
    for nd in path:
        if G.out_degree(nd)==1 and G.in_degree(nd)==1:
            sing_nodes.append(nd)
    return sing_nodes


def get_spades_type_name(count, path, seqs, G, cov=None):
    path_len = len(get_seq_from_path(path,seqs))
    if cov==None:
        cov = get_total_path_mass(path,G)/float(path_len)
    info = ["RNODE", str(count+1), "length", str(path_len),
     "cov", '%.5f' % (cov)]
    return "_".join(info)

def get_contigs_of_mates(node, bamfile, G):
    """ retrieves set of nodes mapped to by read pairs
        having one mate on node; discards isolated nodes
        because they tend to reflect irrelevant alignments
    """
    mate_tigs = set([])
    if node[-1] == "'": node=node[:-1]
    try:    
        for hit in bamfile.fetch(node):
            nref = bamfile.getrname(hit.next_reference_id)
            if nref != node:
                mate_tigs.add(nref)

    except ValueError:
        pass
    source_name = node #re.sub('NODE_','EDGE_', node)

    # print "before removal", mate_tigs
    to_remove = set([])
    for nd in mate_tigs:
        # flip name from "NODE_" prefix back to "EDGE_"
        # differs between contigs set and graph node names
        nd_name = nd #re.sub('NODE_','EDGE_', nd)
        if (G.in_degree(nd_name)==0 and G.out_degree(nd_name)==0) or \
        (not G.has_node(nd_name)):
            to_remove.add(nd)
        # see if nd reachable by node or vice-versa
        # try both flipping to rc and switching source and target    
        elif not any([nx.has_path(G, source_name, nd_name), nx.has_path(G, rc_node(source_name),nd_name), 
          nx.has_path(G, nd_name, source_name), nx.has_path(G, nd_name, rc_node(source_name))]):
            to_remove.add(nd)
    mate_tigs -= to_remove
    # print "after removal", mate_tigs

    return mate_tigs

def is_good_cyc(path, G, bamfile):
    """ check all non-repeat nodes only have mates 
        mapping to contigs in the cycle, ignoring mappings
        to isolated nodes or non-reachable nodes
    """
    sing_nodes = get_non_repeat_nodes(G,path)
    # print sing_nodes
    for nd in sing_nodes:
        mate_tigs = get_contigs_of_mates(nd, bamfile, G)  #re.sub('EDGE_', 'NODE_' ,nd), bamfile, G)
        # print mate_tigs
        # mate_tigs_fixed_names = [re.sub('NODE_','EDGE_', x) for x in mate_tigs]
        # print mate_tigs_fixed_names
        # need to check against F and R versions of path nodes 
        in_path = [x in path for x in mate_tigs] #_fixed_names]
        # print in_path
        path_rc = [rc_node(x) for x in path]
        in_rc_path = [x in path_rc for x in mate_tigs] #_fixed_names]
        # print in_rc_path
        if any([ (not in_path[i] and not in_rc_path[i]) for i in range(len(mate_tigs))]):
            return False
    return True

