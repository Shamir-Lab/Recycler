import subprocess

p = subprocess.Popen("cd ~/BEAR/scripts/",
	shell=True, stdout=subprocess.PIPE)
p.wait()
data_dir = '/home/nasheran/rozovr/recycle_paper_data/'

for size in [400,1600]:#[100,200,400,800,1600]:
	for mult in [2,4]:
		# if size!=100: break
		num_reads = 1250000 * (size/100)
		pref = data_dir+'plasmids_sim_ref_'+str(size)

		####### generate command strings
		abnd_cmmd = 'parametric_abundance.pl ' + pref \
		+'.fasta low > '+pref+'.abnd'
		# print abnd_cmmd 

		gen_reads_cmmd = 'python generate_reads.py -r ' \
		+pref+'.fasta -a '+pref+'.abnd -o ' +pref+\
		'_reads'+str(mult)+'x -t '+str(num_reads)+' -l 100 -i 500 -s 100'
		# print gen_reads_cmmd
		
		zip_cmmd = 'gzip '+ pref+'_reads'+str(mult)+'x.*'
		# print zip_cmmd

		# spades_cmmd = '$eran/SPAdes-3.5.0-Linux/bin/spades.py -1 '\
		# +pref+'_reads.1.fasta.gz -2 '+\
		# pref+'_reads.2.fasta.gz --only-assembler -k 21,33,55 -o '+\
		# data_dir+'ref_'+str(size)+'/'
		# print spades_cmmd

		####### run commands in succession
		# p = subprocess.Popen(abnd_cmmd,
		# shell=True, stdout=subprocess.PIPE)
		# p.wait()
		
		p = subprocess.Popen(gen_reads_cmmd,
		shell=True, stdout=subprocess.PIPE)
		p.wait()
		p = subprocess.Popen(zip_cmmd,
		shell=True, stdout=subprocess.PIPE)
		p.wait()
		