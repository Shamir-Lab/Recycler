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

def get_good_plasmids_dict_from_file(pfile, min_len, max_len, res=None):
    if res==None:
        res = {}
    fp = open(pfile, 'r')
    for name,seq,qual in readfq(fp):
        if min_len <= len(seq) <=max_len and \
        not 'N' in seq:
            res[name] = seq
    return res

########################################################
# Entry point

# read in Aya's plasmid scaffolds, skip any containing N's
# any that are less than 1000 bp
data_dir = '/home/nasheran/rozovr/recycle_paper_data/'
aya_plasmids_file = data_dir + 'Aya_plasmid_scaffolds_PNAS.fa'
ncbi_plasmids_file = data_dir + 'NCBI_plasmids.fasta'

min_len = 1000
max_len = 20000
res = {}
res = get_good_plasmids_dict_from_file(aya_plasmids_file, min_len, max_len, res)
res = get_good_plasmids_dict_from_file(ncbi_plasmids_file, min_len, max_len, res)

for size in [100,200,400,800,1600]:
    samp_file = data_dir + 'plasmids_sim_ref_'+str(size)+'.fasta'
    fo = open(samp_file, 'w')
    sample_inds = random.sample(xrange(len(res)), size)
    keys = res.keys()
    for ind in sample_inds:
    	fo.write(">"+keys[ind]+"\n"+res[keys[ind]]+"\n")
