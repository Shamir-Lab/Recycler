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

def get_cov_from_SPAdes_name(name):
    if 'cov' in G.node[name]:
        return G.node[name]['cov']
    else:
        name_parts = name.split("_")
        cov = name_parts[5]
        if cov[-1]=="'": cov=cov[:-1]
        return float(cov)

def get_SPAdes_base_mass(G, name):
    length = get_length_from_SPAdes_name(name)
    coverage = get_cov_from_SPAdes_name(name)
    return length * coverage

def get_total_path_length(path):
    return sum([get_length_from_SPAdes_name(n) for n in path])

def get_path_coverage_median(path):
    covs = np.array([get_cov_from_SPAdes_name(n) for n in path])
    
    return np.median(covs)

def get_path_coverage_CV(path):
    covs = np.array([get_cov_from_SPAdes_name(n) for n in path])
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

def sort_func(a):
    return a[1]


def rc_node(node):
    """ gets reverse complement
        spades node label
    """ 
    if node[-1] == "'": return node[:-1]
    else: return node + "'"

def get_heaviest_node(nodes, covered_nodes):
    """ given node list sorted by base mass,
        finds first uncovered node
    """
    ind = 0
    while(nodes[ind] in covered_nodes):
        ind+=1
    return nodes[ind]

def enum_high_mass_shortest_paths(G):
    """ given component subgraph, returns list of paths that
        - is non-redundant (includes) no repeats of same cycle
        - includes all cycles starting at node n, taking 
        shortest path (by 1/(length * coverage)) to each of 
        its predecessors, and returning to n
    """
    nodes = []
    nodes[:] = G.nodes()
    unq_sorted_paths = set([])
    paths = []
    # nodes.sort(key=get_SPAdes_base_mass, reverse=True)
    # assign weight edge weights to be 1/mass of starting node
    # making sure node covered whenever rc of node is lets us avoid
    # issue of rc of node having different weight than node
    for e in G.edges():
        G.add_edge(e[0], e[1], cost = 1./get_SPAdes_base_mass(G, e[0]))
        # print e, G[e[0]][e[1]]['cost']

    for node in nodes:
        if node[-1] == "'": continue
        for pred in G.predecessors(node):
            # print pred, get_SPAdes_base_mass(G, pred)
            if pred not in G: continue
            # pred_paths =  [p for p in nx.shortest_path(G, source=highest,
                # target=pred, weight='cost')]
            # for path in pred_paths: # in this case there should only be one path
            try:
                path = nx.shortest_path(G, source=node,
                    target=pred, weight='cost')
            except nx.exception.NetworkXNoPath:
                continue

            all_rc_path = []
            for p in path:
                if p[-1] != "'": p = p+"'"
                all_rc_path.append(p)

            srt = "".join(sorted(all_rc_path))
            if srt not in unq_sorted_paths:
                unq_sorted_paths.add(srt)
                paths.append(tuple(path))
                # for node in path:
                #     covered_nodes.add(node)
                #     covered_nodes.add(rc_node(node))
                # print path
    return paths

def get_path_coverage_bound(path,direction):
    cov = get_cov_from_SPAdes_name(path[0])
    for j in path:
        val = get_cov_from_SPAdes_name(j)
        if direction == 'min' and val < cov:
            cov = val
        elif direction == 'max' and val > cov:
            cov = val
    return cov

def get_total_path_mass(path):
    return sum([get_length_from_SPAdes_name(p) * \
        get_cov_from_SPAdes_name(p) for p in path])

def update_node_coverage_vals(path, comp):
    """ given a path, updates node coverage values
        assuming mean observed path coverage is used
    """
    # print path
    path_copy = list(path)
    tot = get_total_path_mass(path)
    # print "total mass is ", tot
    mean_cov = tot / get_total_path_length(path)
    removed = []
    removed_mean = 0
    # print "mean path coverage is ", mean_cov
    for nd in path:
        nd2 = rc_node(nd)
        
        if comp.in_degree(nd)==1 and comp.out_degree(nd)==1 and\
        comp.in_degree(nd2)==1 and comp.out_degree(nd2)==1:
            
            removed.append(nd)
            # print "removed node ", nd
        
    for nd in removed:
        nd2 = rc_node(nd)
        if nd in comp:
            comp.remove_node(nd)
        if nd2 in comp:
            comp.remove_node(nd2)
        path_copy.remove(nd)
    if len(removed)>0:
        removed_mean = np.mean([get_cov_from_SPAdes_name(p) for p in removed])
    # print [get_cov_from_SPAdes_name(p) for p in removed]
    # print "mean of removed nodes is ", removed_mean
    if removed_mean==0: # no nodes removed - remove all
        new_cov = 1e-6
    elif removed_mean < mean_cov:
        new_cov = mean_cov - removed_mean
    else:
        new_cov = 1e-6

    for nd in path_copy:
        nd2 = rc_node(nd)
        comp.node[nd]['cov'] = new_cov
        if nd2 in comp:
            comp.node[nd2]['cov'] = new_cov

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

def parse_user_input():
    parser = argparse.ArgumentParser(description='recycle extracts cycles likely to be plasmids from metagenome and genome assembly graphs')
    parser.add_argument('-g','--graph', help='(SPAdes 3.50+) FASTG file to process [recommended: before_rr.fastg]',
     required=True, type=str)
    parser.add_argument('-s',
        '--sequences', help='FASTA file (contigs of interest in the graph) to process',
         required=True, type=str)
    parser.add_argument('-l', '--length', help='minimum length required for reporting [default: 1000]',
     required=False, type=int, default=1000)
    parser.add_argument('-m', '--max_CV',
     help='coefficient of variation used for pre-selection [default: 0.50, higher--> less restrictive]; Note: not a requisite for selection',
      required=False, default=1./2, type=float)

    return parser.parse_args()



####### entry point  ##############
# inputs: spades assembly fastg, fasta corresponding
# to component of graph
# outputs: fasta of cycles found by joining component edges

###################################
# 1. read in fastg, create start and end node for each edge
# create edge between start/end, edges out of end using forward 
# strand, edges into start using r.c. strand
args = parse_user_input()
fastg_name = args.graph
comp_name = args.sequences
max_CV = args.max_CV
min_length = args.min_length
fp = open(fastg_name, 'r')
files_dir = os.path.dirname(fp.name)
G = nx.DiGraph() 
lines = get_adj_lines(fastg_name)
G = nx.parse_adjlist(lines, create_using=G)
cycs_ofile = fp.name.replace(".fastg", ".cycs.fasta")
f_cycs_out = open(cycs_ofile, 'w')


###################################
# 2. get subgraph defined by component
# in fasta, remove sources & sinks (can't be in cycle)
fp = open(comp_name, 'r')
comp_nodes = []
seqs = {}
for name,seq,qual in readfq(fp):
    # avoid sources & sinks
    # TODO: get rid of higher order sources/sinks - e.g., 
    # sinks caused by removal of sinks, ...
    if G.out_degree(name)!=0 and G.in_degree(name)!=0:
        comp_nodes.extend([name, name+"'"])
        seqs[name] = seq

if (not (len(comp_nodes)>0 and len(comp_nodes)%2==0)):
    print "graph cleaning implied no cycles possible"
    quit()

comp = G.subgraph(comp_nodes)


###################################
# 3. generate all shortest cycles from each node
# to each of its predecessors until all nodes included
# in some cycle

# naive approach (needs to be optimized for speed):
# remove single node cycs ahead of time

final_paths = set([])
to_remove = set([])
# remove single node cycles, store if long enough

for nd in comp.nodes_with_selfloops(): #nodes_with_selfloops()
    if get_length_from_SPAdes_name(nd) >= min_length:
        if (rc_node(nd),) not in final_paths:
            final_paths.add((nd,))
    to_remove |= set([nd,rc_node(nd)])

for nd in to_remove:
    comp.remove_node(nd)
   

paths = enum_high_mass_shortest_paths(comp)


print "================== path, coverage levels when added ===================="

curr_paths = paths
last_paths = []
while(curr_paths != last_paths):
    last_paths[:] = curr_paths
    curr_paths = []
    if(len(paths)==0): break
    paths.sort(key=get_path_coverage_CV, reverse=False) # low to high
    if get_total_path_mass(paths[0])<1:
        remove_path_nodes_from_graph(paths[0],comp)
        paths.pop()
        curr_paths=paths
        paths = enum_high_mass_shortest_paths(comp)
        continue
    if get_path_coverage_CV(paths[0]) <= max_CV:
        updated_covs = [get_cov_from_SPAdes_name(a) for a in paths[0]]
        curr_paths.append(paths[0])
        covs_before_update = [get_cov_from_SPAdes_name(p) for p in paths[0]]

        update_node_coverage_vals(paths[0], comp)
        # clean_end_nodes_iteratively(comp)
        if get_total_path_length(paths[0])>=min_length:
            final_paths.add(paths[0])
            print paths[0]
            print covs_before_update, "\n"
        paths = enum_high_mass_shortest_paths(comp)

print "==================final_paths identities after updates: ================"
for p in final_paths:
    print p
    print ""

# write out to fasta file
for ind, p in enumerate(final_paths):
    seq = get_seq_from_path(p, seqs)
    info = ["RNODE", str(ind+1), "length", str(len(seq)),
     "cov", '%.5f' % (get_total_path_mass(p)/get_total_path_length(p))]
    # print "_".join(info) + ":", ", ".join(p)
    f_cycs_out.write(">" + "_".join(info) + ": " + ", ".join(p) + "\n" + seq + "\n")



