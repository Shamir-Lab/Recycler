import re, argparse, os
import networkx as nx
import numpy as np
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

def get_length_from_SPAdes_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[3]
    return int(contig_length)

def get_cov_from_SPAdes_name(name,G):
    if name not in G:
        return 0
    if 'cov' in G.node[name]:
        return G.node[name]['cov']
    else:
        name_parts = name.split("_")
        cov = name_parts[5]
        if cov[-1]=="'": cov=cov[:-1]
        return float(cov)

def get_SPAdes_base_mass(G, name):
    length = get_length_from_SPAdes_name(name)
    coverage = get_cov_from_SPAdes_name(name,G)
    return length * coverage



def get_path_coverage_CV(path,G):
    covs = np.array([get_cov_from_SPAdes_name(n,G) for n in path])
    if len(covs)< 2: return 0.000001
    mean = np.mean(covs)
    std = np.std(covs)
    return std/mean


def get_adj_lines(fastg):
    lines = []
    fp = open(fastg, 'r')
    # count = 0
    for name,seq,qual in readfq(fp):
        # count += 1
        # if count % 2 == 0: continue 
        name = re.sub('[:,]'," ", name[:-1])
        lines.append(name)
    return lines

def rc_node(node):
    """ gets reverse complement
        spades node label
    """ 
    if node[-1] == "'": return node[:-1]
    else: return node + "'"

def get_unoriented_sorted_str(path):
    """ creates unq orientation oblivious string rep. of path, 
        used to make sure node covered whenever rc of node is; 
        lets us avoid issue of rc of node having different weight than node
    """
    all_rc_path = []
    for p in path:
        if p[-1] != "'": p = p+"'"
        all_rc_path.append(p)
    return "".join(sorted(all_rc_path))

def enum_high_mass_shortest_paths(G, seen_paths=[]):
    """ given component subgraph, returns list of paths that
        - is non-redundant (includes) no repeats of same cycle
        - includes all cycles starting at node n, taking 
        shortest path (by 1/(length * coverage)) to each of 
        its predecessors, and returning to n
    """

    nodes = []
    nodes[:] = G.nodes()

    unq_sorted_paths = set([])
    # in case orientation obliv. sorted path strings passed in
    for p in seen_paths:
        unq_sorted_paths.add(p)
    paths = []
    # use add_edge to assign edge weights to be 1/mass of starting node
    for e in G.edges():
        G.add_edge(e[0], e[1], cost = 1./get_SPAdes_base_mass(G, e[0]))

    for node in nodes:
        if node[-1] == "'": continue
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


def get_total_path_mass(path,G):
    return sum([get_length_from_SPAdes_name(p) * \
        get_cov_from_SPAdes_name(p,G) for p in path])

def update_node_coverage_vals(path, G, comp, seqs):
    """ given a path, updates node coverage values
        assuming mean observed path coverage is used
    """
    path_copy = list(path)
    tot = get_total_path_mass(path,G)
    mean_cov = tot / get_total_path_length(path, seqs)
    path_determining_nodes = [] # nodes that can be used to est. cov.
    determined_mean = 0

    for nd in path_copy:
        nd2 = rc_node(nd)        
        if G.in_degree(nd)==1 and G.out_degree(nd)==1 and\
            G.in_degree(nd2)==1 and G.out_degree(nd2)==1:
                path_determining_nodes.append(nd)
    if len(path_determining_nodes)>0:
        determined_mean = np.mean([get_cov_from_SPAdes_name(p,G) for p in path_determining_nodes])
        for nd in path_copy:
            nd2 = rc_node(nd)
            if nd in G and nd in comp:
                new_cov = get_cov_from_SPAdes_name(nd,G) - determined_mean
                if new_cov <= 0: 
                    G.remove_node(nd)
                    comp.remove_node(nd)
                else:
                    G.add_node(nd, cov=new_cov)
                    comp.add_node(nd, cov=new_cov)
            if nd2 in G and nd2 in comp:
                new_cov = get_cov_from_SPAdes_name(nd2,G) - determined_mean
                if new_cov <= 0:
                    G.remove_node(nd2)
                    comp.remove_node(nd2)
                else:
                    G.add_node(nd2, cov=new_cov)
                    comp.add_node(nd2, cov=new_cov)
    else: 
        for nd in path_copy:
            nd2 = rc_node(nd)
            if nd in G and nd in comp:
                new_cov = get_cov_from_SPAdes_name(nd,G) - mean_cov
                if new_cov <= 0: 
                    G.remove_node(nd)
                    comp.remove_node(nd)
                else:
                    G.add_node(nd, cov=new_cov)
                    comp.add_node(nd, cov=new_cov)
            if nd2 in G and nd2 in comp:
                new_cov = get_cov_from_SPAdes_name(nd2,G) - mean_cov
                if new_cov <= 0:
                    G.remove_node(nd2)
                    comp.remove_node(nd2)
                else:
                    G.add_node(nd2, cov=new_cov)
                    comp.add_node(nd2, cov=new_cov)
                    

def clean_end_nodes_iteratively(G):
    while(True):
        len_before_update = len(G.nodes())
        for nd in G.nodes():
            nd2 =rc_node(nd)
            if nd not in G or nd not in G: break
            elif G.out_degree(nd)==0 or G.in_degree(nd)==0:
                G.remove_node(nd)
                G.remove_node(nd2)
                break
        if len(G.nodes()) == len_before_update: break

def remove_path_nodes_from_graph(path,G):
    for nd in path:
        if nd in G:
            G.remove_node(nd)
        if rc_node(nd) in G:
            G.remove_node(rc_node(nd))

def get_fasta_stranded_seq(seqs, seq_name):
    """ gets sequence corresponding 
        to same strand as fasta input file 
        or rc seq depending on sequence name
    """
    if seq_name[-1]!="'":
        return seqs[seq_name]
    else: 
        return rc_seq(seqs[seq_name[:-1]])

def get_seq_from_path(path, seqs, max_k_val=55):
    seq = get_fasta_stranded_seq(seqs, path[0])
    if len(path)!=1:
        for p in path[1:]:
            seq += get_fasta_stranded_seq(seqs, p)[max_k_val:]
    return seq

def get_total_path_length(path, seqs):
    # return sum([get_length_from_SPAdes_name(n) for n in path])
    seq = get_seq_from_path(path, seqs)
    return len(seq)

def get_spades_type_name(count, path, seqs, G, cov=None):
    if cov==None:
        cov = get_total_path_mass(path,G)/get_total_path_length(path, seqs)
    info = ["RNODE", str(count+1), "length", str(get_total_path_length(path, seqs)),
     "cov", '%.5f' % (cov)]
    return "_".join(info)

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'recycle extracts cycles likely to be plasmids from metagenome and genome assembly graphs'
        )
    parser.add_argument('-g','--graph',
     help='(SPAdes 3.50+) FASTG file to process [recommended: before_rr.fastg]',
     required=True, type=str
     )
    parser.add_argument('-s',
        '--sequences', help='FASTA file (contigs of interest in the graph) to process',
         required=True, type=str
         )
    parser.add_argument('-l', '--length',
     help='minimum length required for reporting [default: 1000]',
     required=False, type=int, default=1000
     )
    parser.add_argument('-m', '--max_CV',
     help='coefficient of variation used for pre-selection '+
     '[default: 0.50, higher--> less restrictive]; Note: not a requisite for selection',
      required=False, default=1./2, type=float
      )

    return parser.parse_args()



####### entry point  ##############
# inputs: spades assembly fastg, fasta corresponding
# to component of graph
# outputs: fasta of cycles found by joining component edges

###################################
# 1. read in fastg, load graph, create output handle
args = parse_user_input()
fastg_name = args.graph
seqs_name = args.sequences
max_CV = args.max_CV
min_length = args.length
fp = open(fastg_name, 'r')
files_dir = os.path.dirname(fp.name)
G = nx.DiGraph() 
lines = get_adj_lines(fastg_name)
G = nx.parse_adjlist(lines, create_using=G)
# fasta of sequences
fasta_ofile = fp.name.replace(".fastg", ".cycs.fasta")
f_cycs_fasta = open(fasta_ofile, 'w')
# file containing path name (corr. to fasta), path, coverage levels
# when path is added
cycs_ofile = fp.name.replace(".fastg", ".cycs.paths_w_cov.txt")
f_cyc_paths = open(cycs_ofile, 'w')

###################################
# 2. get subgraph defined by component fasta
# remove sources & sinks (can't be in cycle)
fp = open(seqs_name, 'r')
seq_nodes = []
seqs = {}
for name,seq,qual in readfq(fp):
    # avoid sources & sinks
    # TODO: get rid of higher order sources/sinks - e.g., 
    # sinks caused by removal of sinks, ...
    if G.out_degree(name)!=0 and G.in_degree(name)!=0:
        seq_nodes.extend([name, name+"'"])
        seqs[name] = seq

if (not (len(seq_nodes)>0 and len(seq_nodes)%2==0)):
    print "graph cleaning implied no cycles possible"
    quit()

G = G.subgraph(seq_nodes)


###################################
# 3. generate all shortest cycles from each node
# to each of its predecessors 

# remove single node cycs ahead of time

self_loops = set([])
non_self_loops = set([])
final_paths_dict = {}
to_remove = set([])
# remove single node cycles, store if long enough

path_count = 0
for nd in G.nodes_with_selfloops(): #nodes_with_selfloops()
    if get_length_from_SPAdes_name(nd) >= min_length:
        if (rc_node(nd),) not in self_loops:
            name = get_spades_type_name(path_count, (nd,), seqs, G)
            self_loops.add((nd,))
            final_paths_dict[name] = (nd,)
            path_count += 1
    to_remove |= set([nd,rc_node(nd)])

for nd in to_remove:
    if nd in G: G.remove_node(nd)

# H = G.to_undirected()

print "================== path, coverage levels when added ===================="
for comp in list(nx.strongly_connected_component_subgraphs(G)):
    # comp = G.subgraph(comp.nodes())
 
    # initialize shortest path set considered
    paths = enum_high_mass_shortest_paths(comp)

    # peeling - iterate until no change in path set from 
    # one iteration to next

    
    last_path_count = 0
    last_node_count = 0
    # continue as long as you either removed a low mass path
    # from the component or added a new path to final paths
    while(path_count!=last_path_count or\
        len(comp.nodes())!=last_node_count):
        last_node_count = len(comp.nodes())
        last_path_count = path_count

        if(len(paths)==0): break

        # print "paths ", paths 
        # using initial set of paths or set from last iteration
        # sort the paths by CV and test the lowest CV path for removal 
        # need to use lambda because get_cov needs 2 parameters
        sort_fun = lambda x: get_cov_from_SPAdes_name(x,G)
        paths.sort(key=sort_fun, reverse=False) # low to high
        curr_path = paths[0]
        if get_total_path_mass(curr_path,G)<1:
            remove_path_nodes_from_graph(curr_path,comp)
            # recalc. paths since some might have been disrupted
            non_self_loops.add(get_unoriented_sorted_str(curr_path))
            paths = enum_high_mass_shortest_paths(comp,non_self_loops)
            continue
        if get_path_coverage_CV(curr_path,G) <= max_CV and \
        get_unoriented_sorted_str(curr_path) not in non_self_loops:
            
            covs_before_update = [get_cov_from_SPAdes_name(p,G) for p in curr_path]
            cov_val_before_update = get_total_path_mass(curr_path,G) /\
             get_total_path_length(curr_path, seqs)
            update_node_coverage_vals(curr_path, G, comp, seqs)
            # clean_end_nodes_iteratively(comp)
            name = get_spades_type_name(path_count,curr_path, seqs, G, cov_val_before_update)
            path_count += 1
            non_self_loops.add(get_unoriented_sorted_str(curr_path))
            final_paths_dict[name] = curr_path

            # only report to file if long enough
            if get_total_path_length(curr_path, seqs)>=min_length:
                print curr_path
                print "before", covs_before_update
                print "after", [get_cov_from_SPAdes_name(p,G) for p in curr_path]
                print len(comp.nodes()), " nodes remain in component\n"
                f_cyc_paths.write(name + "\n" +str(curr_path)+ "\n" + str(covs_before_update) + "\n")
            # recalculate paths on the component
            paths = enum_high_mass_shortest_paths(comp,non_self_loops)

# done peeling
# print final paths to screen
print "==================final_paths identities after updates: ================"

# write out sequences to fasta
for p in final_paths_dict.keys():
    seq = get_seq_from_path(final_paths_dict[p], seqs)
    print final_paths_dict[p]
    print ""
    if len(seq)>=min_length:
        f_cycs_fasta.write(">" + p + "\n" + seq + "\n")    

