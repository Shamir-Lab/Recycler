#!/usr/bin/env python

import argparse, os
from recyclelib.utils import *
import pysam

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'Recycler extracts likely plasmids (and other circular DNA elements) from de novo assembly graphs'
        )
    parser.add_argument('-g','--graph',
     help='(spades 3.50+) assembly graph FASTG file to process; recommended for spades 3.5: before_rr.fastg, for spades 3.6+:assembly_graph.fastg',
     required=True, type=str
     )
    parser.add_argument('-k','--max_k',
        help='integer reflecting maximum k value used by the assembler',
        required=True, type=int, default=55
        )
    parser.add_argument('-b','--bam',
        help='BAM file resulting from aligning reads to contigs file, filtering for best matches',
        required=True, type=str
        )
    parser.add_argument('-l', '--length',
     help='minimum length required for reporting [default: 1000]',
     required=False, type=int, default=1000
     )
    parser.add_argument('-m', '--max_CV',
     help='coefficient of variation used for pre-selection [default: 0.5, higher--> less restrictive]',
      required=False, default=1./2, type=float
      )
    parser.add_argument('-i','--iso',
        help='True or False value reflecting whether data sequenced was an isolated strain',
        required=False, type=bool, default=False
        )
    parser.add_argument('-o','--output_dir',
        help='Output directory',
        required=False, type=str
        )
    return parser.parse_args()


if __name__ == '__main__':

    ####### entry point  ##############
    # inputs: spades assembly fastg, BAM of reads aligned to unitigs
    # outputs: fasta of cycles found by joining component edges

    ###################################
    # read in fastg, load graph, create output handle
    args = parse_user_input()
    fastg = args.graph
    max_CV = args.max_CV
    max_k = args.max_k
    min_length = args.length
    fp = open(fastg, 'r')
    files_dir = os.path.dirname(fp.name)

    # output 1 - fasta of sequences
    if args.output_dir:
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        basename = os.path.basename(fp.name)
        out_fp = os.path.join(args.output_dir, basename)
        (root, ext) = os.path.splitext(out_fp)
    else:
        (root,ext) = os.path.splitext(fp.name)

    fasta_ofile = root + ext.replace(".fastg", ".cycs.fasta")
    f_cycs_fasta = open(fasta_ofile, 'w')
    # output 2 - file containing path name (corr. to fasta),
    # path, coverage levels when path is added
    cycs_ofile = root + ext.replace(".fastg", ".cycs.paths_w_cov.txt")
    f_cyc_paths = open(cycs_ofile, 'w')
    bamfile = pysam.AlignmentFile(args.bam)
    ISO = args.iso

    ###################################
    # graph processing begins


    G = get_fastg_digraph(fastg)
    G.remove_nodes_from(list(nx.isolates(G)))

    cov_vals = [get_cov_from_spades_name(n) for n in G.nodes()]
    MED_COV = np.median(cov_vals)
    STD_COV = np.std(cov_vals)
    # set thresholds for max. CV, min
    # path coverage for allowing cross mappings
    if ISO:
        thresh = np.percentile(cov_vals, 95)
    else:
        thresh = np.percentile(cov_vals, 75)

    print(MED_COV, STD_COV, thresh)
    path_count = 0
    # gets set of long simple loops, removes short
    # simple loops from graph
    SEQS = get_fastg_seqs_dict(fastg,G)
    long_self_loops = get_long_self_loops(G, min_length, SEQS)
    non_self_loops = set([])
    VISITED_NODES = set([]) # used to avoid problems due to RC nodes we may have removed
    final_paths_dict = {}

    for nd in long_self_loops:
        name = get_spades_type_name(path_count, nd,
            SEQS, G, get_cov_from_spades_name(nd[0]))
        final_paths_dict[name] = nd
        path_count += 1

    comps = nx.strongly_connected_component_subgraphs(G)
    COMP = nx.DiGraph()
    redundant = False
    print("================== path, coverage levels when added ====================")

    ###################################
    # iterate through SCCs looking for cycles
    for c in comps:
        # check if any nodes in comp in visited nodes
        # if so continue
        for node in c.nodes():
            if c in VISITED_NODES:
                redundant = True
                break
        if redundant:
            redundant = False
            continue # have seen the RC version of component
        COMP = c.copy()

        # initialize shortest path set considered
        paths = enum_high_mass_shortest_paths(COMP)

        # peeling - iterate until no change in path set from
        # one iteration to next

        last_path_count = 0
        last_node_count = 0
        # continue as long as you either removed a low mass path
        # from the component or added a new path to final paths
        while(path_count!=last_path_count or\
            len(COMP.nodes())!=last_node_count):

            last_node_count = len(COMP.nodes())
            last_path_count = path_count

            if(len(paths)==0): break

            # using initial set of paths or set from last iteration
            # sort the paths by CV and test the lowest CV path for removal
            # need to use lambda because get_cov needs 2 parameters

            # make tuples of (CV, path)
            path_tuples = []
            for p in paths:
                path_tuples.append((get_wgtd_path_coverage_CV(p,COMP,SEQS,max_k_val=max_k), p))

            # sort in ascending CV order
            path_tuples.sort(key=lambda path: path[0])

            curr_path = path_tuples[0][1]
            # print paths
            # print curr_path
            if get_unoriented_sorted_str(curr_path) not in non_self_loops:
                path_mean, _ = get_path_mean_std(curr_path, G, SEQS, max_k_val=max_k)

                ## only report to file if long enough and good
                ## first good case - paired end reads on non-repeat nodes map on cycle
                ## typical or low coverage level
                ## second good case - high coverage (pairs may map outside due to high chimericism),
                ## near constant coverage level
                if (
                    len(get_seq_from_path(curr_path, SEQS, max_k_val=max_k))>=min_length \
                    and is_good_cyc(curr_path,G,bamfile) and \
                    get_wgtd_path_coverage_CV(curr_path,COMP,SEQS,max_k_val=max_k) <= (max_CV/len(curr_path))
                    ) or \
                (
                    len(get_seq_from_path(curr_path, SEQS, max_k_val=max_k))>=min_length and (path_mean > thresh) \
                    and get_wgtd_path_coverage_CV(curr_path,COMP,SEQS,max_k_val=max_k) <= (max_CV/len(curr_path))
                    ):
                    print(curr_path)
                    non_self_loops.add(get_unoriented_sorted_str(curr_path))
                    name = get_spades_type_name(path_count, curr_path, SEQS, COMP)
                    # covs = [get_cov_from_spades_name_and_graph(p,COMP) for p in curr_path]
                    covs = get_path_covs(curr_path,COMP)
                    print("before", covs)
                    f_cyc_paths.write(name + "\n" +str(curr_path)+ "\n" + str(covs)
                        + "\n" + str([get_num_from_spades_name(p) for p in curr_path]) + "\n")
                    update_path_coverage_vals(curr_path, COMP, SEQS)
                    path_count += 1
                    print("after", get_path_covs(curr_path,COMP))
                    final_paths_dict[name] = curr_path

                # recalculate paths on the component
                print(len(COMP.nodes()), " nodes remain in component\n")

                paths = enum_high_mass_shortest_paths(COMP,non_self_loops)
        rc_nodes = [rc_node(n) for n in COMP.nodes()]
        VISITED_NODES.update(COMP.nodes())
        VISITED_NODES.update(rc_nodes)

    # done peeling
    # print final paths to screen
    print("==================final_paths identities after updates: ================")

    # write out sequences to fasta
    for p in final_paths_dict.keys():
        seq = get_seq_from_path(final_paths_dict[p], SEQS, max_k_val=max_k)
        print(final_paths_dict[p])
        print(" ")
        if len(seq)>=min_length:
            f_cycs_fasta.write(">" + p + "\n" + seq + "\n")
