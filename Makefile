SPADES_DIR = /home/nasheran/rozovr/SPAdes-3.6.2-Linux/bin
BWA_PATH = /usr/bin/bwa
SAMTOOLS_PATH = /usr/bin/samtools
# paired end fastq files here
READS_DIR = /home/nasheran/rozovr/recycle_paper_data
# spades outputs including assembly graph here, output gets written here
NUM_THREADS = 32
REF_CNT = 100
ASSEM_PREF = ref_$(REF_CNT)
ASSEMBLY_DIR = $(READS_DIR)/$(ASSEM_PREF)
READ_F1 = plasmids_sim_ref_$(REF_CNT)_reads.1.fasta.gz
READ_F2 = plasmids_sim_ref_$(REF_CNT)_reads.2.fasta.gz

tests:
	nosetests --nocapture

all: assemble index_and_map filter_to_primary sort_and_index recycle clean

assemble:
	# $(SPADES_DIR)/spades.py -t $(NUM_THREADS) --only-assembler -k 21,33,55 -1 $(READS_DIR)/$(READ_F1) -2 $(READS_DIR)/$(READ_F2) -o $(ASSEMBLY_DIR) 
	$(SPADES_DIR)/spades.py -t $(NUM_THREADS) -1 $(READS_DIR)/$(READ_F1) -2 $(READS_DIR)/$(READ_F2) -o $(ASSEMBLY_DIR) 
	python make_fasta_from_fastg.py -g $(ASSEMBLY_DIR)/assembly_graph.fastg

index_and_map: $(ASSEMBLY_DIR)/assembly_graph.nodes.fasta
	$(BWA_PATH) index $^
	$(BWA_PATH) mem -t $(NUM_THREADS) $^ $(READS_DIR)/$(READ_F1) $(READS_DIR)/$(READ_F2) | $(SAMTOOLS_PATH) view -buS - > $(READS_DIR)/$(ASSEM_PREF)_pe.bam

filter_to_primary: $(READS_DIR)/$(ASSEM_PREF)_pe.bam
	$(SAMTOOLS_PATH) view -bF 0x0800 $^  > $(READS_DIR)/$(ASSEM_PREF)_pe_primary.bam

sort_and_index: $(READS_DIR)/$(ASSEM_PREF)_pe_primary.bam
	$(SAMTOOLS_PATH) sort $^ $(READS_DIR)/$(ASSEM_PREF)_pe_primary.sort
	$(SAMTOOLS_PATH) index $(READS_DIR)/$(ASSEM_PREF)_pe_primary.sort.bam

recycle: $(READS_DIR)/$(ASSEM_PREF)_pe_primary.sort.bam
	python recycle.py -g $(ASSEMBLY_DIR)/assembly_graph.fastg -b $^

clean:
	rm $(READS_DIR)/$(ASSEM_PREF)_pe.bam $(READS_DIR)/$(ASSEM_PREF)_pe_primary.bam