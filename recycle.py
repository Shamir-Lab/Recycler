import argparse, os
from recycle.utils import *
import pysam


def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'recycle extracts cycles likely to be plasmids from metagenome and genome assembly graphs'
        )
    parser.add_argument('-g','--graph',
     help='(spades 3.50+) FASTG file to process [recommended: before_rr.fastg]',
     required=True, type=str
     )
    parser.add_argument('-l', '--length',
     help='minimum length required for reporting [default: 1000]',
     required=False, type=int, default=1000
     )
    parser.add_argument('-m', '--max_CV',
     help='coefficient of variation used for pre-selection '+
     '[default: 0.25, higher--> less restrictive]; Note: not a requisite for selection',
      required=False, default=1./4, type=float
      )
    # parser.add_argument('-b','--bamp', 
    #     help='prefix to BAM files resulting from aligning reads to fasta file', 
    #     required=False, type=str
    #     )
    
    return parser.parse_args()



####### entry point  ##############
# inputs: spades assembly fastg, BAM of reads aligned to unitigs
# outputs: fasta of cycles found by joining component edges

###################################
# 1. read in fastg, load graph, create output handle
args = parse_user_input()
fastg = args.graph
max_CV = args.max_CV
min_length = args.length
fp = open(fastg_name, 'r')
files_dir = os.path.dirname(fp.name)

# output 1 - fasta of sequences
(root,ext) = os.path.splitext(fp.name)
fasta_ofile = root + ext.replace(".fastg", ".cycs.fasta")
f_cycs_fasta = open(fasta_ofile, 'w')
# output 2 - file containing path name (corr. to fasta), 
# path, coverage levels when path is added
cycs_ofile = root + ext.replace(".fastg", ".cycs.paths_w_cov.txt")
f_cyc_paths = open(cycs_ofile, 'w')


###################################
# 2a. extract self loop edges from nodes having
# AND-types == outies on both ends at most 500 bp away from end

# ands_file = args.bamp + '.fasta.ands.srt.bam'
# samfile = pysam.AlignmentFile(ands_file)

# print "before adding, ", len(G.edges()), " edges"

# for node in G.nodes():
#     try:
#         hits = samfile.fetch(node)
#         num_hits = sum(1 for _ in hits)

#         # print len(hits), hits 
#         if num_hits>1:
#             G.add_edge(node,node)
#             G.add_edge(rc_node(node),rc_node(node))
#     except ValueError:
#         continue

# print "after adding, ", len(G.edges()), " edges"


# next use contig_joining_type to connect
# sink nodes having more than one pair of reads
# connecting them
# joins_file = args.bamp + '.fasta.joins.srt.bam'
# samfile = pysam.AlignmentFile(joins_file)
# print "before adding, ", len(G.edges()), " edges"

# sinks = []
# hit_cnts = {}
# for node in G.nodes():
#     if G.out_degree(node)==0:
#         sinks.append(node)
# # print len(sinks), " sinks: ", sinks
# for node in sinks:
#     try:
#         hits = samfile.fetch(node)
        
#         for hit in samfile.fetch(node):
#             nref = samfile.getrname(hit.next_reference_id)
            
#             if nref in sinks: 
#                 hit_cnts[(node,nref)] = hit_cnts.get((node,nref),0)+1
#                 if hit_cnts[(node,nref)]>1:
#                     G.add_edge(node, nref)
#                     G.add_edge(rc_node(nref), rc_node(node))

#     except ValueError:
#         continue
            
# print "after adding, ", len(G.edges()), " edges"


# 2b. get subgraph defined by component fasta
# remove sources & sinks (can't be in cycle)
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

# if (not (len(seq_nodes)>0 and len(seq_nodes)%2==0)):
#     print "graph cleaning implied no cycles possible"
#     quit()

# G = G.subgraph(seq_nodes)


###################################
# 3. generate all shortest cycles from each node
# to each of its predecessors 


G = get_fastg_digraph(fastg)

# gets set of long simple loops, removes short
# simple loops from graph
long_self_loops = get_long_self_loops(G, min_length)
VISITED_NODES = set([]) # used to avoid problems due to RC nodes we may have removed

comps = nx.strongly_connected_component_subgraphs(G)
COMP = nx.DiGraph()

print "================== path, coverage levels when added ===================="
for c in comps:
    # check if any nodes in comp in visited nodes
    # if so continue
    
    COMP = c.copy()
    SEQS = get_fastg_seqs_dict(fastg, COMP)

    # initialize shortest path set considered
    paths = enum_high_mass_shortest_paths(COMP)

    
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

        # print paths
        # using initial set of paths or set from last iteration
        # sort the paths by CV and test the lowest CV path for removal 
        # need to use lambda because get_cov needs 2 parameters
        
        # make tuples of (CV, path)
        path_tuples = []
        for p in paths:
            # path_tuples.append((get_path_coverage_CV(p,G), p))
            # print get_wgtd_path_coverage_CV(p,G,seqs), p
            path_tuples.append((get_wgtd_path_coverage_CV(p,G,seqs), p))
        
        # sort in ascending CV order
        path_tuples.sort(key=lambda path: path[0]) 
        
        curr_path = path_tuples[0][1]
        if get_total_path_mass(curr_path,G)<1:
            
            remove_path_nodes_from_graph(curr_path,comp)
            # recalc. paths since some might have been disrupted
            non_self_loops.add(get_unoriented_sorted_str(curr_path))
            paths = enum_high_mass_shortest_paths(comp,non_self_loops)
            continue

        # if get_path_coverage_CV(curr_path,G) <= max_CV and \
        if get_wgtd_path_coverage_CV(curr_path,G,seqs) <= max_CV and \
        get_unoriented_sorted_str(curr_path) not in non_self_loops:

            covs_before_update = [get_cov_from_spades_name_and_graph(p,G) for p in curr_path]
            cov_val_before_update = get_total_path_mass(curr_path,G) /\
             get_total_path_length(curr_path, seqs)
            path_nums = [get_num_from_spades_name(p) for p in curr_path]
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
                print "after", [get_cov_from_spades_name_and_graph(p,G) for p in curr_path]
                f_cyc_paths.write(name + "\n" +str(curr_path)+ "\n" + str(covs_before_update) 
                    + "\n" + str(path_nums) + "\n")
            # recalculate paths on the component
            print len(comp.nodes()), " nodes remain in component\n"
            
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

