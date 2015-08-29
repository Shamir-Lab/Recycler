all: fetch_scripts

clean:
	rm fetch_*

fetch_scripts: fetch_joins fetch_ands

fetch_joins: extract_contig_joining_pairs.c
	gcc -o $@ $^ -L. -lbam -lz

fetch_ands: extract_mate_pair_type_reads.c
	gcc -o $@ $^ -L. -lbam -lz