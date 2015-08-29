#include <stdio.h>
#include "bam.h"
#include "sam.h"

// to compile, in ~/plasmid_code dir run: gcc -o fetch_joins extract_contig_joining_type_reads.c -L. -lbam -lz
// based on https://www.biostars.org/p/77802/

int main(int argc, char* argv[]) {
    samfile_t *ifile = NULL, *ofile = NULL;
    bam1_t *read = bam_init1();
    int keep = 0;
    char *p = NULL;

    //Open input file, either SAM or BAM
    p = strrchr(argv[1], '.');
    if(strcmp(p, ".bam") == 0) {
        ifile = samopen(argv[1], "rb", NULL);
    } else {
        ifile = samopen(argv[1], "r", NULL);
    }

    bam_header_t *head = ifile->header;

    //Open output file
    // ofile = samopen("contig_joining_type.bam", "wb", ifile->header);
    ofile = samopen(argv[2], "wb", ifile->header);

    //Iterate through the lines
    while(samread(ifile, read) > 1) {
        keep = 0;
        //Is the read's mate on a different chromosome/contig?
        if(read->core.tid != read->core.mtid) {
            // are both mates mapped?
            if(!(read->core.flag & BAM_FUNMAP) && !(read->core.flag & BAM_FMUNMAP)){ 
               keep = 1;
            }
            
        }
        if(keep) samwrite(ofile, read);
    }
    bam_destroy1(read);
    samclose(ifile);
    samclose(ofile);
    return 0;
}