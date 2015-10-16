SCRIPTS_DIR = ~/plasmid_scripts
DATA_DIR = /home/nasheran/rozovr/recycle_paper_data
REF_CNT = 100
INPUT_DIR = $(DATA_DIR)/ref_$(REF_CNT)
INPUT = $(INPUT_DIR)/before_rr

all: $(INPUT).1000.fasta


# filter contigs by length
$(INPUT).1000.fasta: $(INPUT).fasta
	perl $(SCRIPTS_DIR)/remove.smalls.pl 1000 $^ > $@

# I. using perl scripts, find contigs having matching ends
# 1) split contigs in middle, results put in new directory chopContig
# .PHONY: chopped
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
$(INPUT_DIR)/circular.fna: $(INPUT_DIR)/fasta
	perl $(SCRIPTS_DIR)/collection.pl $(INPUT_DIR)/fasta > $@




# # II. using perl scripts, find contigs having reads mapping to opposite ends
# # 1) Extract both ends (500 bp) of the original contig (more than 1kb)
# perl $(SCRIPTS_DIR)/fasta_manipulate.Endseq.pl ../ayabrown/trimmed_reads/U_res/contigs4.1000.fasta 500
# # It will generate 2 files: contig.leftseq.fa and contig.rightseq.fa
# # 2) make fasta from reads fastqs
# awk 'NR % 4 == 1 {print ">" $0 } NR % 4 == 2 {print $0}' ../ayabrown/trimmed_reads/U_1_trimmed.fastq > ../ayabrown/trimmed_reads/U_1_trimmed.fasta & 
# awk 'NR % 4 == 1 {print ">" $0 } NR % 4 == 2 {print $0}' ../ayabrown/trimmed_reads/U_2_trimmed.fastq > ../ayabrown/trimmed_reads/U_2_trimmed.fasta &
# # 3) The name of raw Illumina PE reads need to be modified: remove blank and add length values to the ends of the head line
# perl $(SCRIPTS_DIR)/fasta_manipulate.PEheader.pl ../ayabrown/trimmed_reads/U_1_trimmed.fasta > ../ayabrown/trimmed_reads/U_1_trimmed.name.fasta &
# perl $(SCRIPTS_DIR)/fasta_manipulate.PEheader.pl ../ayabrown/trimmed_reads/U_2_trimmed.fasta > ../ayabrown/trimmed_reads/U_2_trimmed.name.fasta &
# # 4) Preparation of 2 databases for blastn
# formatdb -i contig.leftseq.fa -p F -o T
# formatdb -i contig.rightseq.fa -p F -o T
# # 5) BLASTing (takes ~25 mins) per process, 
# blastall -p blastn -d contig.leftseq.fa -i ../ayabrown/trimmed_reads/U_1_trimmed.name.fasta -o read2contigLeft.blastn -e 1e-10 -F F -m 8
# blastall -p blastn -d contig.rightseq.fa -i ../ayabrown/trimmed_reads/U_2_trimmed.name.fasta -o read2contigRight.blastn -e 1e-10 -F F -m 8
# # [To define the number of Illumina short reads (100 bp, 150 bp, 250 bp) matched to each end of each contigs, we used the following criteria: 
# # 100% alignment identity, alignment length of hit read at least 90 bp, and alignment length accounting for the hit read sequence should be greater than 99%.]
# perl $(SCRIPTS_DIR)/parse_PEread.1.pl read2contigRight.blastn > hitread2Right.tab
# perl $(SCRIPTS_DIR)/parse_PEread.1.pl read2contigLeft.blastn > hitread2Left.tab
# # For the left end of the contig, only reserve the reads that are complementary reverse to the contig sequence.
# perl $(SCRIPTS_DIR)/parse_PEread.2.pl hitread2Left.tab > hitread2Left.reverse.tab
# # If one of a read pair located on the left of the contig, check the other one of a read pair whether located on the right of the contig or not. 
# # If true, it will result in 2 files: the contig name (circular.PE.list) and the corresponding read names (circular.PE.read.list).
# perl $(SCRIPTS_DIR)/parse_PEread.3.pl hitread2Left.reverse.tab hitread2Right.tab
# # 6) Remove the duplicated contig names and reserve only one per contig
# sort circular.PE.list | uniq > circular.PE.sort.list
# # Use the intersection of 2 tests as final result
# perl $(SCRIPTS_DIR)/intersection.pl circular.fna circular.PE.sort.list > circular.final.list
# # extract fasta file of results
# $(SCRIPTS_DIR)/faSomeRecords circular.fna circular.final.list circular.final.fa
