#!/usr/bin/env python

## Script generates an alignment of variable sites, a list of positions
## Allows for masking of given elements of genome

from Bio import SeqIO
import os, gzip, sys

from optparse import OptionParser

usage = '''
cleanAlignment instructions
Usage: cleanAlignment.py [options] alignment_folder

alignment_folder: folder containing output from getAlignment


Options
-v, --varsites  proportion of variable sites to be called for variable site to be retained
				default 0.70
-s, --seq 		proportion of sequence (variable sites) required to be called for a 
				sequence to be retained, default 0.70
-n, --align		length of alignment before start filtering out sequences, default 0

Example:

./cleanAlignment.py -v 0.70 -s 0.70 alignment_folder


'''


if __name__=='__main__':

	### PROCESS INPUT OPTIONS AND ARGUMENTS ###	
	parser = OptionParser()
	parser.add_option( '-v', '--varsites', action = 'store', type='float', dest = 'varsite_keep', default = 0.70 )
	parser.add_option( '-s', '--seq', action = 'store', type='float', dest = 'seq_keep', default = 0.70 )
	parser.add_option( '-n', '--align', action = 'store', type='int', dest = 'align_n', default = 0 )
	opts, args = parser.parse_args()
	
	#check have correct number of arguments
	if len(args) == 1:
		prefix = args[0]
		varsite_keep = opts.varsite_keep
		seq_keep = opts.seq_keep
		align_n = opts.align_n
	else:
		sys.stdout.write(usage)
		sys.exit(1)
	
	
	sys.stdout.write('Cleaning %s\n'%prefix)

	fafile = '%s_snps.fa'%prefix
	posfile = '%s_positions.txt'%prefix
	rejectfile = '%s_reject.txt'%prefix

	f = open(posfile, 'r')
	pos = [l.strip() for l in f]
	f.close()

	bases = 'ACGT'

	seq = SeqIO.parse(fafile, 'fasta')
	seqlist = [s for s in seq]
	n = len(seqlist)

	if n==1:
		clean_aln = '%s_clean_snps.fa'%prefix
		SeqIO.write( seqlist , clean_aln, "fasta" )
		clean_posfile = '%s`'%prefix
		f = open(clean_posfile, 'w')
		f.close()
		sys.stdout.write('Done\n\n')

	if n>1:
	
		#get list of positions called >= varsite_keep of time
		seq_generator = zip(*seqlist)
		keep = [i for ( i, a ) in enumerate( seq_generator ) 
					   if float(len([ai for ai in a if ai in bases]))/n >=varsite_keep]

		#write clean posfile
		clean_posfile = '%s_clean_positions.txt'%prefix
		f = open(clean_posfile, 'w')
		for i,p in enumerate(pos):
			if i in keep:
				f.write('%s\n'%p)

		f.close()

		#clean alignment
		for seq in seqlist:
			clean_bases = "".join( seq.seq[ i ] for i in keep )
			seq.seq._data = clean_bases

		#search new alignment for any sequence not called at >=seq_keep of sites
		n = len(keep)
		reject = []
		rejectDict = {}
		for s in seqlist:
			called = len([b for b in s.seq if b in bases])
			if n>0 and n>align_n and float(called)/n < seq_keep:
				reject.append(s.id)
				rejectDict[s.id] = [n, called, float(called)/n]



		if len(reject)>0:
			sys.stdout.write('\tReject: %s\n'%", ".join(reject))
			rf = open(rejectfile, 'w')
			for r in list(rejectDict.keys()):
				rf.write('%s\t%s\t%s\t%0.4f\n'%(r, rejectDict[r][0], rejectDict[r][1], rejectDict[r][2]))
			rf.close()
		
		
		#write intermediate alignment with excess sites removed, but all samples in still
		int_aln = '%s_pre_exclude_snps.fa'%prefix
		seqlist = [seq for seq in seqlist]
		SeqIO.write( seqlist , int_aln, "fasta" )
		
		
		#clean alignment - write file with samples excluded
		clean_aln = '%s_clean_snps.fa'%prefix
		seqlist = [seq for seq in seqlist if seq.id not in reject]
		SeqIO.write( seqlist , clean_aln, "fasta" )

		sys.stdout.write('Done\n\n')
		sys.stdout.flush()