# creates a report on the number of plasmids reported
# and their alignment relative to some reference 
INPUT_DIR = /vol/scratch/rozovr/M_res
REPEAT_RES_1 = 0
REPEAT_RES_2 = 0
REF = $(INPUT_DIR)/previously_validated_M_contigs.fa

ifeq ($(REPEAT_RES_1),0)
	INPUT1 = before_rr.fasta
else
	INPUT1 = contigs.fasta
endif

ifeq ($(REPEAT_RES_2),0)
	INPUT2 = before_rr
else
	INPUT2 = contigs
endif


all: $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).cycs.dbl.fasta $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.summary

# 1) concat output cycles to self:
$(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).cycs.dbl.fasta: $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).cycs.fasta
	python ~/recycle/concat_seq_to_self.py -i $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).cycs.fasta

# 2) run nucmer
$(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer: $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).cycs.dbl.fasta
	echo "great"
	/home/gaga/rozovr/MUMmer3.23/nucmer $(REF) $^ -p $@

# 3) parse alignments, write out summary 
$(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.summary: $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.delta 
	touch $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.delta | \
	awk '$10==100.00 && $15==100.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.delta | \
	awk '$10==100.00 && $15>=90.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.summary
	/home/gaga/rozovr/MUMmer3.23/show-coords -r -c -l $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.delta | \
	awk '$10==100.00 && $15>=80.00' | cut -d'|' --complement -f 1-6 | uniq | wc -l >> $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.summary

clean:
	rm $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).cycs.dbl.fasta
	rm $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.delta
	rm $(INPUT_DIR)/iter2_on_$(INPUT1)/$(INPUT2).nucmer.summary