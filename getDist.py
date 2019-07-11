#!/usr/bin/env python


import sys, re, os, gzip, time, multiprocessing, math, random
from Bio import SeqIO

from optparse import OptionParser

usage = '''
getDist instructions
Usage: getDist.py [options] alignment.fa output_prefix

alignment.fa: multifasta file alignment of variable sites, e.g. from getAlignment.py
output_prefix: output folder and filename stem

Options
-p, --nprocs  number of cores to use, default 1


Example:

./getDist.py -p 15 alignment.fa output_prefix


'''

bases = "ACGT"

def get_distance( seq1, seq2 ):
	'''Function to calculate the [Hamming] distance between two sequences'''
	return sum(c1!=c2 for c1, c2 in zip( seq1, seq2 ) if c1 in bases and c2 in bases)

if __name__=="__main__":
	
	### PROCESS INPUT OPTIONS AND ARGUMENTS ###	
	parser = OptionParser()
	parser.add_option( '-p', '--nprocs', action = 'store', type='int', dest = 'nprocs', default = 1 )
	opts, args = parser.parse_args()
	
	#check have correct number of arguments
	if len(args) == 2:
		inname, outname = args
		nprocs = opts.nprocs
	else:
		sys.stdout.write(usage)
		sys.exit(1)
	
	## open multifasta of variable sites
	seqlist = list(SeqIO.parse(inname, "fasta"))
	sys.stdout.write("Successfully opened snps multi-fasta file.\n")
	sys.stdout.flush()
	
	##shuffle seqlist
	random.shuffle(seqlist)
	
	## Do the data matrix
	start = time.perf_counter()
	
	def compare (i, seqlist, start_i, end_i, out_q):
		'''compare function compares a subset of the data'''
		outfile = "%s_%s"%(outname, i)
		# print([s.id for s in seqlist])  # debug
		f = open(outfile, "w")
		for i in seqlist[start_i:end_i]:
			for j in seqlist:
				if i.id < j.id:
					d = get_distance(i.seq, j.seq)
					f.write("%s\t%s\t%s\n"%(i.id, j.id, d))
	
	
	
	out_q = multiprocessing.Queue()
	chunksize = int(math.ceil(len(seqlist) / float(nprocs)))
	procs = []
	
	for i in range(nprocs):
		p = multiprocessing.Process(
				target=compare,
				args=(i, seqlist, i*chunksize, (i+1)*chunksize, out_q))
		procs.append(p)
		p.start()
	
	
	# Wait for all worker processes to finish
	for p in procs:
		p.join()
	
	
	time_elapsed = time.perf_counter() - start
	sys.stdout.write("Successfully wrote distance matrix in %s seconds.\n"%time_elapsed)
	sys.stdout.write("Done.\n\n")
