#include <stdio.h>
#include "bam.h"
#include "sam.h"

// to compile, in ~/plasmid_code dir run: gcc -o filter_ors extract_outie_ends.c -L. -lbam -lz
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
    ofile = samopen("OR_type_filtered.bam", "wb", ifile->header);

    //Iterate through the lines
    while(samread(ifile, read) > 1) {
        keep = 0;
        //Is the read's mate on the same chromosome/contig?
        if(read->core.tid == read->core.mtid) {
            //Are the mates on opposite strands?
            if(read->core.flag & BAM_FREVERSE && !(read->core.flag & BAM_FMREVERSE)) {
                if(read->core.pos < read->core.mpos) {
                    // Is one of the mates 500 bp or less from an end?
                    if (read-> core.pos <= 500 || read->core.mpos > head->target_len[read->core.tid] - 500)
                        keep=1;
                }
            } else if(!(read->core.flag & BAM_FREVERSE) && read->core.flag & BAM_FMREVERSE) {
                if(read->core.mpos < read->core.pos) {
                    if (read-> core.mpos <= 500 || read->core.pos > head->target_len[read->core.tid] - 500)
                        keep=1;
                }
            }
        }
        if(keep) samwrite(ofile, read);
    }
    bam_destroy1(read);
    samclose(ifile);
    samclose(ofile);
    return 0;
}