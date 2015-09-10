
INPUT_DIR = /vol/scratch/rozovr/M_res/
READS_DIR = /home/gaga/rozovr/recycle_paper_data/


all: fetch_scripts


#### fetch scripts used for processing bam file ####
# output bams used either for 2 step assembly or input
# to recycler
fetch_scripts: fetch_joins fetch_ands

fetch_joins: extract_contig_joining_pairs.c
	gcc -o $@ $^ -L. -lbam -lz

fetch_ands: extract_mate_pair_type_reads.c
	gcc -o $@ $^ -L. -lbam -lz


















clean:
	rm fetch_*
