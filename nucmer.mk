# creates a report on the number of plasmids reported
# and their alignment relative to some reference 
var = 5

ifeq ($(REPEAT_RES_1),0)
	INPUT1 = before_rr.fasta
else
	INPUT1 = contigs.fasta
endif

ifeq ($(REPEAT_RES_2),0)
	INPUT2 = before_rr.fasta
	GRAPH2 = before_rr.fastg
else
	INPUT2 = contigs.fasta
	GRAPH2 = contigs.fastg
endif

all:
	echo "the var is $(var)"

# 1) concat output cycles to self:
# python ~/recycle/concat_seq_to_self.py -i $eran/recycle_paper_data/ref_800/before_rr.cycs.fasta
# [will also need to remove illegal chars from reference, as done in simulation]

# 2) run nucmer
# $eran/MUMmer3.23/nucmer $eran/recycle_paper_data/plasmids_sim_ref_100.cln.fasta $eran/recycle_paper_data/ref_100/before_rr.cycs.dbl.fasta -p $eran/recycle_paper_data/ref_100/nucmer                                              

# 3) parse alignments for 100/100, making sure aligned pairs are unique
# $eran/MUMmer3.23/show-coords -r -c -l $eran/recycle_paper_data/ref_100/nucmer.delta | awk '$10==100.00 && $15==100.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l
# $eran/MUMmer3.23/show-coords -r -c -l $eran/recycle_paper_data/ref_100/nucmer.delta | awk '$10==100.00 && $15>=90.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l
# $eran/MUMmer3.23/show-coords -r -c -l $eran/recycle_paper_data/ref_100/nucmer.delta | awk '$10==100.00 && $15>=80.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l
