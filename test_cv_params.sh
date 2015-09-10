#!/bin/sh
set -e # exit on first error
for ref in 100 200 400 800 1600; do
	for cv in 0.25 0.375 0.5; do
		make -f nucmer.sim.mk clean CV="$cv" REF_CNT="$ref"
		make -f nucmer.sim.mk all CV="$cv" REF_CNT="$ref"
	done
done

# def update_node_coverage_vals(path, G, comp, seqs, max_k_val=55):
#     """ given a path, updates node coverage values
#         assuming mean observed path coverage is used
#     """
#     path_copy = list(path)
#     # tot = get_total_path_mass(path,G)
#     # mean_cov = tot / get_total_path_length(path, seqs)
#     # path_determining_nodes = [] # nodes that can be used to est. cov.
#     # determined_mean = 0

#     # for nd in path_copy:
#     #     nd2 = rc_node(nd)        
#     #     if G.in_degree(nd)==1 and G.out_degree(nd)==1 and\
#     #         G.in_degree(nd2)==1 and G.out_degree(nd2)==1:
#     #             path_determining_nodes.append(nd)
#     # if len(path_determining_nodes)>0:
#     #     determined_mean = np.mean([get_cov_from_SPAdes_name(p,G) for p in path_determining_nodes])
#     #     for nd in path_copy:
#     #         nd2 = rc_node(nd)
#     #         if nd in G and nd in comp:
#     #             new_cov = get_cov_from_SPAdes_name(nd,G) - determined_mean
#     #             if new_cov <= 0: 
#     #                 G.remove_node(nd)
#     #                 comp.remove_node(nd)
#     #             else:
#     #                 G.add_node(nd, cov=new_cov)
#     #                 comp.add_node(nd, cov=new_cov)
#     #         if nd2 in G and nd2 in comp:
#     #             new_cov = get_cov_from_SPAdes_name(nd2,G) - determined_mean
#     #             if new_cov <= 0:
#     #                 G.remove_node(nd2)
#     #                 comp.remove_node(nd2)
#     #             else:
#     #                 G.add_node(nd2, cov=new_cov)
#     #                 comp.add_node(nd2, cov=new_cov)
#     # else: 
#     covs = np.array([get_cov_from_SPAdes_name(n,G) for n in path])
#     # if len(covs)< 2: return 0.000001
#     # mean = np.mean(covs)
#     wgts = np.array([(get_length_from_SPAdes_name(n)-max_k_val) for n in path])
#     tot_len = get_total_path_length(path, seqs)
#     wgts = np.multiply(wgts, 1./tot_len)
#     mean_cov = np.average(covs, weights = wgts)

#     for nd in path_copy:
#         nd2 = rc_node(nd)
#         if nd in G and nd in comp:
#             new_cov = get_cov_from_SPAdes_name(nd,G) - mean_cov
#             if new_cov <= 0: 
#                 G.remove_node(nd)
#                 comp.remove_node(nd)
#             else:
#                 G.add_node(nd, cov=new_cov)
#                 comp.add_node(nd, cov=new_cov)
#         if nd2 in G and nd2 in comp:
#             new_cov = get_cov_from_SPAdes_name(nd2,G) - mean_cov
#             if new_cov <= 0:
#                 G.remove_node(nd2)
#                 comp.remove_node(nd2)
#             else:
#                 G.add_node(nd2, cov=new_cov)
#                 comp.add_node(nd2, cov=new_cov)
#     # clean_end_nodes_iteratively(G)
#     # clean_end_nodes_iteratively(comp)    
