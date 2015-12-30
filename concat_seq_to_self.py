import argparse, os

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

def parse_user_input():
    parser = argparse.ArgumentParser(
        description=
        'given multi-fasta of contigs, creates new file with each sequence concatenated to itself'
        )
    parser.add_argument('-i','--input',
     help='multi FASTA file to process ', required=True, type=str
     )
    parser.add_argument('-m','--maxk',
     help='maximum value of k used in assembly, needed to avoid repeats in glueing ', required=False, type=int, default=55
     )

    return parser.parse_args()

args = parse_user_input()
fasta_name = args.input


fp = open(fasta_name, 'r')
(root,ext) = os.path.splitext(fp.name)
out_name = root + ext.replace(".fasta", ".dbl.fasta")
fo = open(out_name, 'w')

for name,seq,qual in readfq(fp):
    dbl_seq = seq + seq[args.maxk:]
    fo.write(">"+name+'\n'+dbl_seq+"\n")
