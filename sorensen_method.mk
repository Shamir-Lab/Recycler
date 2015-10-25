SCRIPTS_DIR = ~/plasmid_scripts
DATA_DIR = /home/nasheran/rozovr/recycle_paper_data
REF_CNT = 100
INPUT_DIR = $(DATA_DIR)/ref_$(REF_CNT)
INPUT = $(INPUT_DIR)/before_rr
REF = $(INPUT_DIR)/../plasmids_sim_ref_$(REF_CNT).cln.fasta



all: $(INPUT_DIR)/circular.1000.fasta chopped $(INPUT_DIR)/minimus2.out $(INPUT_DIR)/circular.fasta
validate: $(INPUT_DIR)/circular.dbl.fasta $(INPUT_DIR)/circular.nucmer $(INPUT_DIR)/circular.nucmer.summary

# filter contigs by length
$(INPUT).1000.fasta: $(INPUT).fasta
	perl $(SCRIPTS_DIR)/remove.smalls.pl 1000 $^ > $@

# I. using perl scripts, find contigs having matching ends
# 1) split contigs in middle, results put in new directory chopContig
chopped: $(INPUT).1000.fasta
	perl $(SCRIPTS_DIR)/chop.sequence.pl $^

# 2) Merge overlapping sequences using minimus2 with default parameters
# 3) Put all the resulting *.fasta files into a new folder fasta. 
# Overlapping sequences at the previous ends have been removed.
$(INPUT_DIR)/minimus2.out: $(INPUT_DIR)/chopContig
	perl $(SCRIPTS_DIR)/minimus2.pl $^ > $@
	mkdir $(INPUT_DIR)/fasta
	cp $(INPUT_DIR)/chopContig/*.fasta $(INPUT_DIR)/fasta/

# 4) Put all detected circular contigs into a file
$(INPUT_DIR)/circular.fasta: $(INPUT_DIR)/fasta
	perl $(SCRIPTS_DIR)/collection.pl $(INPUT_DIR)/fasta > $@

# ------------------- validation steps -----------------------

$(INPUT_DIR)/circular.dbl.fasta: $(INPUT_DIR)/circular.fasta
	python ~/recycle/concat_seq_to_self.py -m 0 -i $^

$(INPUT_DIR)/circular.nucmer: $(INPUT_DIR)/circular.dbl.fasta
	/home/gaga/rozovr/MUMmer3.23/nucmer $(REF) $^ -p $@

# 3) parse alignments, write out summary 
$(INPUT_DIR)/circular.nucmer.summary: $(INPUT_DIR)/circular.nucmer.delta 
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT_DIR)/circular.nucmer.delta | \
	awk '$$10==100.00 && $$15==100.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l > $(INPUT_DIR)/circular.nucmer.summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT_DIR)/circular.nucmer.delta | \
	awk '$$10==100.00 && $$15>=90.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT_DIR)/circular.nucmer.summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT_DIR)/circular.nucmer.delta | \
	awk '$$10==100.00 && $$15>=80.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT_DIR)/circular.nucmer.summary

