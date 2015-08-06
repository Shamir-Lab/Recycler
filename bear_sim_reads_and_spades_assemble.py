bear_dir =  '~/BEAR' 
data_dir = '/home/nasheran/rozovr/recycle_paper_data/'

for size in [100,200,400,800,1600]:
	num_reads = 1250000 * (size/100)
	pref = data_dir+'plasmids_sim_ref_'+str(size)
	abnd_cmmd = 'parametric_abundance.pl ' + pref \
	+'.fasta low > '+pref+'.abnd'
	print abnd_cmmd 

	gen_reads_cmmd = 'python generate_reads.py -r ' \
	+pref+'.fasta -a '+pref+'.abnd -o ' +pref+\
	'_reads -t '+str(num_reads)+' -l 100 -i 500 -s 100'
	print gen_reads_cmmd
	
	zip_cmmd = 'gzip '+ pref+'_reads.*'
	print zip_cmmd

	spades_cmmd = '$eran/SPAdes-3.5.0-Linux/bin/spades.py -1 '\
	+pref+'_reads.1.fasta.gz -2 '+\
	pref+'_reads.2.fasta.gz --only-assembler -k 21,33,55 -o '+\
	data_dir+'ref_'+str(size)+'/'
	print spades_cmmd

