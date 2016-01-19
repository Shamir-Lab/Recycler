# gets simple (single contig) cycles from 
# (usually plasmid) metagenomes, leaves rest of 
# graph as is

import re, os, argparse

def readfq(fp): # this is a generator function
    """ # lh3's fast fastX reader: 
        https://github.com/lh3/readfq/blob/master/readfq.py
    """
    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def get_length_from_SPAdes_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[3]
    return int(contig_length)

def parse_user_input():
    parser = argparse.ArgumentParser(description='gets simple (single contig) cycles from plasmid metagenomes, leaves rest of graph as is; outputs these to two separate files')
    parser.add_argument('-i','--input', help='Input (SPAdes 3.50+) FASTG to process', required=True, type=str)
    parser.add_argument('-m','--min_length', help='Minimum cycle length to keep (shorter cycles put in new graph file; default = 1000)', 
        required=False, type=int, default=1000)
    return parser.parse_args()

###########################################
# ENTRY POINT
# inputs: fastg, min length (def 1000)
# outputs: fasta containing simple cycs, 
# fastg containing rest of graph
args = parse_user_input()
fastg_name = args.input
min_length = args.min_length

cycs = set([])
fp = open(fastg_name, 'r')

simple_cycs_ofile = fp.name.replace(".fastg", ".simple.cycs.fasta")
f_cycs_out = open(simple_cycs_ofile, 'w')
filt_graph_ofile = fp.name.replace(".fastg", ".simple.cycs_filtered.fastg")
f_graph_out = open(filt_graph_ofile, 'w')


#### potential issue: output fastg 
#### may contain links to nodes that have been removed
#### may be potential hazard downstream
cnt = 0
for name,seq,qual in readfq(fp):
    cnt+=1
    if cnt%2!=1:
        continue
    orig_name = name
    name = re.sub('[:,]'," ", name[:-1])
    name_parts = name.split(" ")
    
    if (len(name_parts)==2) and (name_parts[0]==name_parts[1]):
        length = get_length_from_SPAdes_name(name_parts[0])
        if length >= min_length and name_parts[0][:-1] not in cycs:
            # second condition to avoid forward & rc being included together
            # non-apostrophe version always comes first in fastg
            cycs.add(name_parts[0])
            f_cycs_out.write(">"+name_parts[0]+"\n"+seq+"\n")
    elif name_parts[0] not in cycs and name_parts[0][:-1] not in cycs:
        f_graph_out.write(">"+orig_name+"\n"+seq+"\n")


f_cycs_out.close()
f_graph_out.close()

