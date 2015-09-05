import numpy as np
complements = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}


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

def rc_seq(dna):
    rev = reversed(dna)
    return "".join([complements[i] for i in rev])

def get_num_from_SPAdes_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[1]
    return int(contig_length)      

def get_length_from_SPAdes_name(name):
    name_parts = name.split("_")
    contig_length = name_parts[3]
    return int(contig_length)

def get_cov_from_SPAdes_name(name,G):
    if name not in G:
        return 0
    if 'cov' in G.node[name]:
        return G.node[name]['cov']
    else:
        name_parts = name.split("_")
        cov = name_parts[5]
        if cov[-1]=="'": cov=cov[:-1]
        return float(cov)

def get_SPAdes_base_mass(G, name):
    length = get_length_from_SPAdes_name(name)
    coverage = get_cov_from_SPAdes_name(name,G)
    return length * coverage

def get_seq_from_path(path, seqs, max_k_val=55):
    seq = get_fasta_stranded_seq(seqs, path[0])
    if len(path)!=1:
        for p in path[1:]:
            seq += get_fasta_stranded_seq(seqs, p)[max_k_val:]
    return seq

def get_total_path_length(path, seqs):
    # return sum([get_length_from_SPAdes_name(n) for n in path])
    seq = get_seq_from_path(path, seqs)
    return len(seq)

def get_wgtd_path_coverage_CV(path, G, seqs, max_k_val=55):
    covs = np.array([get_cov_from_SPAdes_name(n,G) for n in path])
    if len(covs)< 2: return 0.000001
    # mean = np.mean(covs)
    wgts = np.array([(get_length_from_SPAdes_name(n)-max_k_val) for n in path])
    mean = np.average(covs, weights = wgts)
    tot_len = get_total_path_length(path, seqs)
    wgts = np.multiply(wgts, 1./tot_len)   
    std = sum(np.dot(wgts,(covs-mean)**2))

    # if mean == 0: return 1000 
    return std/mean


def get_path_coverage_CV(path,G):
    covs = np.array([get_cov_from_SPAdes_name(n,G) for n in path])
    if len(covs)< 2: return 0.000001
    mean = np.mean(covs)
    std = np.std(covs)
    # if mean == 0: return 1000 
    return std/mean

def get_total_path_mass(path,G):
    return sum([get_length_from_SPAdes_name(p) * \
        get_cov_from_SPAdes_name(p,G) for p in path])

def get_fasta_stranded_seq(seqs, seq_name):
    """ gets sequence corresponding 
        to same strand as fasta input file 
        or rc seq depending on sequence name
    """
    if seq_name[-1]!="'":
        return seqs[seq_name]
    else: 
        return rc_seq(seqs[seq_name[:-1]])