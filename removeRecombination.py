#!/usr/bin/env python

#############################################################################
# removeRecombination.py
# remove recombination from an alignment using CFML output
# 28 Mar 2017
#############################################################################


from optparse import OptionParser
from Bio import SeqIO
import os, sys, re
import subprocess, shlex

usage = '''
removeRecombination instructions
Usage: removeRecombination.py [options] cluster_folder

Example:

bin/removeRecombination.py -c 1 /well/bag/deyre/data/cdif/he 

Options
-c, --cluster		cluster number

'''


if __name__=="__main__":
	
	
	### PROCESS INPUT OPTIONS AND ARGUMENTS ###
	parser = OptionParser()
	parser.add_option( "-c", "--cluster", action = "store", type='int', dest = "cluster", default = 0 )	
	opts, args = parser.parse_args()
	
	#check have correct number of arguments
	if len(args) == 1:
		stem = args[0]
	else:
		sys.stdout.write(usage)
		sys.exit(1)
	
	#save options	
	c = opts.cluster

	#make output directory
	if not os.path.exists('%s/recomb_corr'%stem):
		os.mkdir('%s/recomb_corr'%stem)
	
	#check if CFML output for this cluster
	if not os.path.exists('%s/cluster_ml/cluster_%s_cf.importation_status.txt'%(stem, c)):
		sys.stderr.write('No files found matching input values, please check and try again\n')
		sys.exit()
	
	sys.stdout.write('Removing recombination from cluster %s\n'%c)
	sys.stdout.write('\t- reading inputs\n')
	sys.stdout.flush()
	
	#read in variable sites
	posfile = '%s/cluster/cluster_%s_clean_positions.txt'%(stem, c)
	snpsfile = '%s/cluster/cluster_%s_clean_snps.fa'%(stem, c)
	
	f = open(posfile, 'r')
	pos_all = [int(l.strip()) for l in f]
	f.close()
	
	## if no positions, i.e. no differences between the sequences then exit
	if len(pos_all) == 0:
		sys.stdout.write('\t- no variable sites, exiting\n')
		sys.exit()
	
	seqlist = [s for s in SeqIO.parse(snpsfile, 'fasta')]
	
	## if less than min_cluster_size --> exit
	min_cluster_size = 3
	if len(seqlist) < min_cluster_size:
		sys.stdout.write('\t- not enough samples in this cluster, exiting\n')
		sys.exit()
	
	#read in CFML recombination sites
	impfile = '%s/cluster_ml/cluster_%s_cf.importation_status.txt'%(stem, c)
	f = open(impfile, 'r')
	recomb = [l.strip().split() for l in f][1:]
	recomb = [[int(r[1]), int(r[2])] for r in recomb]
	f.close()
	
	sys.stdout.write('\t- defining positions to keep\n')
	#strip recombining sites from alignment
	remove_pos = []
	for i, p in enumerate(pos_all):
		for r in recomb:
			if p <= r[1] and p >= r[0]:
				#remove position
				remove_pos.append(i)
				break
	
	keep_pos = [i for i in range(len(pos_all)) if i not in remove_pos]
	
	sys.stdout.write('\t- updating sequence\n')
	for s in seqlist:
		s.seq._data = ''.join( s.seq[ i ] for i in keep_pos )
		s.id = s.id.split(' ')[0]
		s.description = ''
	
	sys.stdout.write('\t- writing clean alignment\n')
	#write recombination corrected alignment of variable sites
	SeqIO.write(seqlist, '%s/recomb_corr/cluster_%s_norecomb_snps.fa'%(stem,c), 'fasta')
	#write recombination correct positions
	f = open('%s/recomb_corr/cluster_%s_norecomb_positions.txt'%(stem,c), 'w')
	f.write('\n'.join( str(pos_all[ i ]) for i in keep_pos ))
	f.close()
