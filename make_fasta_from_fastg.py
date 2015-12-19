import re
from utils import *

def get_adj_lines(fastg):
    lines = []
    fp = open(fastg, 'r')
    count = 0
    for name,seq,qual in readfq(fp):
        count += 1
        if count % 2 == 0: continue 
        name = re.sub('[:,]'," ", name[:-1]).split(" ")[0]
        line = ">"+name+"\n"+seq
        print line
    #     lines.append(name)
    # return lines

fastg_name = '/home/nasheran/rozovr/megahit_test/Sheerli_samp/k99.fastg'
# fg = open(fastg_name, 'r')
get_adj_lines(fastg_name)