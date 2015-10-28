# creates a report on the number of plasmids reported
# and their alignment relative to some reference 
INPUT_DIR = /home/nasheran/rozovr/aya/spades3.5_test/M_res/
TWO_STEP_ASSEM = 0
REPEAT_RES_1 = 0
REPEAT_RES_2 = 0
REF = $(INPUT_DIR)/previously_validated_M_contigs.fa

ifeq ($(REPEAT_RES_1),0)
	INPUT1 = before_rr
else
	INPUT1 = contigs
endif

ifeq ($(REPEAT_RES_2),0)
	INPUT2 = before_rr
else
	INPUT2 = contigs
endif

ifeq ($(TWO_STEP_ASSEM),1)
	INPUT = $(INPUT_DIR)/iter2_on_$(INPUT1).fasta/$(INPUT2)
else
	INPUT = $(INPUT_DIR)/$(INPUT1)
endif


all: $(INPUT).cycs.dbl.fasta $(INPUT).nucmer $(INPUT).nucmer.summary

# 1) concat output cycles to self:
$(INPUT).cycs.dbl.fasta: $(INPUT).cycs.fasta
	python ~/recycle/concat_seq_to_self.py -i $(INPUT).cycs.fasta

# 2) run nucmer
$(INPUT).nucmer: $(INPUT).cycs.dbl.fasta
	echo "great"
	/home/gaga/rozovr/MUMmer3.23/nucmer $(REF) $^ -p $@

# 3) parse alignments, write out summary 
$(INPUT).nucmer.summary: $(INPUT).nucmer.delta 
	touch $(INPUT).nucmer.summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT).nucmer.delta | \
	awk '$$10==100.00 && $$15==100.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT).nucmer.summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT).nucmer.delta | \
	awk '$$10==100.00 && $$15>=90.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT).nucmer.summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT).nucmer.delta | \
	awk '$$10==100.00 && $$15>=80.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT).nucmer.summary

clean:
	rm $(INPUT).cycs.dbl.fasta
	rm $(INPUT).nucmer.delta
	rm $(INPUT).nucmer.summary