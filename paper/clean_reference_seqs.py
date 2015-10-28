import argparse
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
        "removes characters not in 'ACGT' from sequences"
        )
    parser.add_argument('-i','--input',
     help='multi-FASTA reference file to process',
     required=True, type=str
     )
    

    return parser.parse_args()

args = parse_user_input()
fasta_name = args.input
fp = open(fasta_name, 'r')
out_name = fp.name.replace(".fasta", ".cln.fasta")
fo = open(out_name, 'w')

for name,seq,qual in readfq(fp):
    clean_seq = []
    for c in seq:
        if c in 'ACGT':
            clean_seq.append(c)
    clean_seq =''.join(clean_seq)
    fo.write('>'+name+"\n"+clean_seq+"\n")

