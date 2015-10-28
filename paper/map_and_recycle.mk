
# INPUT_DIR = /home/nasheran/ayabrown/trimmed_reads/M_res/

READS_DIR = /home/nasheran/rozovr/recycle_paper_data/
REF_CNT = 100
INPUT_DIR = $(READS_DIR)/ref_$(REF_CNT)
INPUT1 = before_rr.fasta

all: fetch_joins fetch_ands $(INPUT_DIR)/$(INPUT1).bwt \
	$(INPUT_DIR)/$(INPUT1).ann $(INPUT_DIR)/$(INPUT1).amb \
	$(INPUT_DIR)/$(INPUT1).sa $(INPUT_DIR)/$(INPUT1) \
	$(INPUT_DIR)/reads_to_$(INPUT1).bam \
	$(INPUT_DIR)/reads_to_$(INPUT1).joins.bam \
	$(INPUT_DIR)/reads_to_$(INPUT1).joins.srt.bam \
	$(INPUT_DIR)/reads_to_$(INPUT1).ands.srt.bam \
	$(INPUT_DIR)/reads_to_$(INPUT1).ands.bam


#### fetch scripts used for processing bam file ####
# output bams used either for 2 step assembly or input
# to recycler

fetch_joins: extract_contig_joining_pairs.c
	gcc -o $@ $^ -L. -lbam -lz

fetch_ands: extract_mate_pair_type_reads.c
	gcc -o $@ $^ -L. -lbam -lz

$(INPUT_DIR)/$(INPUT1).bwt $(INPUT_DIR)/$(INPUT1).ann $(INPUT_DIR)/$(INPUT1).amb $(INPUT_DIR)/$(INPUT1).sa: $(INPUT_DIR)/$(INPUT1)
	~/bwa/bwa index $^

$(INPUT_DIR)/reads_to_$(INPUT1).bam: $(INPUT_DIR)/$(INPUT1).bwt $(INPUT_DIR)/$(INPUT1).ann $(INPUT_DIR)/$(INPUT1).amb $(INPUT_DIR)/$(INPUT1).sa
	~/bwa/bwa mem -t 16 $(INPUT_DIR)/$(INPUT1) \
	$(READS_DIR)/plasmids_sim_ref_$(REF_CNT)_reads.1.fasta.gz \
	$(READS_DIR)/plasmids_sim_ref_$(REF_CNT)_reads.2.fasta.gz | \
	samtools view -buS - > $@

$(INPUT_DIR)/reads_to_$(INPUT1).joins.bam: $(INPUT_DIR)/reads_to_$(INPUT1).bam
	fetch_joins $^ $@

$(INPUT_DIR)/reads_to_$(INPUT1).joins.srt.bam: $(INPUT_DIR)/reads_to_$(INPUT1).joins.bam
	samtools sort $(INPUT_DIR)/reads_to_$(INPUT1).joins.bam $(INPUT_DIR)/reads_to_$(INPUT1).joins.srt
	samtools index $@

$(INPUT_DIR)/reads_to_$(INPUT1).ands.bam: $(INPUT_DIR)/reads_to_$(INPUT1).bam
	fetch_ands $^ $@

$(INPUT_DIR)/reads_to_$(INPUT1).ands.srt.bam: $(INPUT_DIR)/reads_to_$(INPUT1).ands.bam
	samtools sort $(INPUT_DIR)/reads_to_$(INPUT1).ands.bam $(INPUT_DIR)/reads_to_$(INPUT1).ands.srt
	samtools index $@

clean:
	rm -f $(INPUT_DIR)/reads_to_$(INPUT1)*.bam
	rm -f fetch_joins fetch_ands

