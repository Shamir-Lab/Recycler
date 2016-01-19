# creates a report on the number of plasmids reported
# and their alignment relative to some reference
REF_CNT = 100
CV = 0.5
PARENT_DIR = /home/nasheran/rozovr/recycle_paper_data
INPUT_DIR = $(PARENT_DIR)/ref_$(REF_CNT)
# TWO_STEP_ASSEM = 0
# REPEAT_RES_1 = 0
# REPEAT_RES_2 = 0
REF = $(PARENT_DIR)/plasmids_sim_ref_$(REF_CNT).cln.fasta
INPUT = $(INPUT_DIR)/contigs
# ifeq ($(REPEAT_RES_1),0)
# 	INPUT1 = before_rr
# else
# 	INPUT1 = contigs
# endif

# ifeq ($(REPEAT_RES_2),0)
# 	INPUT2 = before_rr
# else
# 	INPUT2 = contigs
# endif

# ifeq ($(TWO_STEP_ASSEM),1)
# 	INPUT = $(INPUT_DIR)/iter2_on_$(INPUT1).fasta/$(INPUT2)
# else
# 	INPUT = $(INPUT_DIR)/$(INPUT1)
# endif

# $(INPUT).cycs.fasta
all:  $(INPUT).cycs.dbl.fasta $(INPUT).nucmer.delta $(INPUT).nucmer$(CV).summary $(INPUT).$(CV).tally

# 0) run recycler on input
# $(INPUT).cycs.fasta: $(INPUT).fastg
# 	python ~/recycle/recycle.py -g $(INPUT).fastg -s $(INPUT).fasta -m $(CV) #-b $(INPUT_DIR)/reads_to_$(INPUT1)

# 1) concat output cycles to self:
$(INPUT).cycs.dbl.fasta: $(INPUT).cycs.fasta
	python ~/recycle/paper/concat_seq_to_self.py -i $(INPUT).cycs.fasta -m 55

# 2) run nucmer
$(INPUT).nucmer.delta: $(INPUT).cycs.dbl.fasta
	/home/gaga/rozovr/MUMmer3.23/nucmer $(REF) $^ -p $(INPUT).nucmer

# 3) parse alignments, write out summaries
$(INPUT).nucmer$(CV).summary $(INPUT).hit_nodes$(CV).txt:  $(INPUT).nucmer.delta 
	grep '>' -c $(INPUT).cycs.fasta > $(INPUT).nucmer$(CV).summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT).nucmer.delta | \
	awk '$$10==100.00 && $$15==100.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT).nucmer$(CV).summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT).nucmer.delta | \
	awk '$$10==100.00 && $$15>=90.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT).nucmer$(CV).summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT).nucmer.delta | \
	awk '$$10==100.00 && $$15>=80.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT).nucmer$(CV).summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT).nucmer.delta  | \
	awk '$$10==100.00 && $$15>=80.00' | cut -d'|' --complement -f 1-6 | cut -f2 > $(INPUT).hit_nodes$(CV).txt

# 4) tally hits by different ranges/properties
$(INPUT).$(CV).tally: $(INPUT).cycs.fasta $(INPUT).cycs.paths_w_cov.txt $(INPUT).hit_nodes$(CV).txt
	python ~/recycle/paper/tally_hits.py -p $(INPUT) -n $(INPUT).hit_nodes$(CV).txt > $@

clean:
	# rm -f $(INPUT).cycs.fasta
	rm -f $(INPUT).cycs.dbl.fasta
	rm -f $(INPUT).nucmer.delta
	rm -f $(INPUT).nucmer$(CV).summary
	rm -f $(INPUT).hit_nodes$(CV).txt
	rm -f $(INPUT).$(CV).tally