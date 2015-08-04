import random
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

########################################################
# Entry point

plasmids_file = '/home/nasheran/rozovr/recycle_paper_data/NCBI_plasmids.fasta'
samp_file = 'samp_plasmids_out.fasta'
max_len = 20000
samp_size = 800
res = {}
print plasmids_file
fp = open(plasmids_file, 'r')
for name,seq,qual in readfq(fp):
	if len(seq)<max_len:
		res[name] = seq

# fo = open(samp_file, 'w')
sample_inds = random.sample(xrange(len(res)), samp_size)
keys = res.keys()
for ind in sample_inds:
	print ">"+keys[ind]+"\n"+res[keys[ind]]+"\n"
	# fo.write(">"+keys[ind]+"\n"+res[keys[ind]]+"\n")
