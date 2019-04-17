#!/usr/bin/env pypy

## Script generates an alignment of variable sites, a list of positions
## Allows for masking of given elements of genome

## uses pypy rather than python 

import sys, re, os, gzip, time, multiprocessing, math, random, datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

from optparse import OptionParser
from subprocess import call

usage = '''
mtAlign instructions
Usage: mtAlign.py [options] list_of_fasta.txt ref.fa output_prefix

list_of_fasta.txt: tab delimited 1) name, 2) fasta path and filename
ref.fa: reference in fasta format
output_prefix: output folder

Options
-m, --mask  optional mask file for mobile elements, 
            tab delimited - start and end coordinates, counting from one, not zero
-p, --nprocs  number of cores to use, default 1

Example:

./mtAlign.py -m mask.txt -p 4 list_of_fasta.txt ref.fa output_prefix


'''



if __name__=='__main__':

    start = datetime.datetime.now()

    ### PROCESS INPUT OPTIONS AND ARGUMENTS ### 
    parser = OptionParser()
    parser.add_option( '-m', '--mask', action = 'store', type='string', dest = 'maskfile', default = '' )
    parser.add_option( '-p', '--nprocs', action = 'store', type='int', dest = 'nprocs', default = 1 )
    opts, args = parser.parse_args()
    
    #check have correct number of arguments
    if len(args) == 3:
        listoffasta, refpath, outname_prefix = args
        maskfile = opts.maskfile
        nprocs = opts.nprocs
    else:
        sys.stdout.write(usage)
        sys.exit(1)
    
    #read in reference, concatenate multiple chromasomes
    #refseq = SeqIO.read( refpath, 'fasta' ).seq
    refseq = "".join([s.seq._data for s in SeqIO.parse( refpath, 'fasta' )])
    reflen = len( refseq )
    sys.stderr.write('Successfully read in reference, length %s.\n'%reflen)    
    
    #read in mask file if provided
    # set up default which is to use all sites, i.e. mask none
    masksites = [1]*reflen
    if maskfile:
        r = open(maskfile,'r')
        mask = [(int(l.split('\n')[0].split('\t')[0])-1, int(l.split('\n')[0].split('\t')[1])-1) for l in r] #convert back to base0
        r.close()
    
        for m in mask:
            for i in range(m[0], m[1]+1):
                masksites[i]=0
    
        sys.stderr.write('Successfully read in mobile elements mask.\n') 
    
    
    ## Read in all the sequences, and replace the id with the required nicename
    seqlist = []
    with open( listoffasta ) as fp:
        for line in fp:
            if line:
                nicename, fapath = line.strip().split('\t')
                if os.path.exists(fapath):
                    sys.stdout.write('%s\n'%fapath)
                    with gzip.open(fapath, 'rt') as fa_fh:
                        fa = list(SeqIO.parse(fa_fh, 'fasta'))[0]
                        fa.id = nicename
                        fa.description = ''
                    seqlist.append(fa)
                else:
                    sys.stderr.write(fapath+' does not exist, skipping...')
    
    time_elapsed = datetime.datetime.now() - start
    time_elapsed = float(time_elapsed.seconds) + float(time_elapsed.microseconds/1e6)
    sys.stdout.write("Successfully read in sequences in %0.1f seconds.\n"%time_elapsed)
    
    
    ## obtain non-shared positions
    start = datetime.datetime.now()
    
    def subset_pos(j, seqlist_j, masksites_j, start_i, end_i, out_q, outname_prefix):
        sys.stdout.write('Starting from %s to %s\n'%(start_i, end_i-1))
        sys.stdout.flush()
        masksites_j = masksites_j[start_i: end_i]
        seqlist_j = [s[start_i: end_i] for s in seqlist_j]
        seq_generator = zip( masksites_j, *seqlist_j  )
        nonshared_pos =[ i+start_i for ( i, a ) in enumerate( seq_generator ) 
                   if a[0]==1 and len( set( [ ai for ai in a[1:] if ai in 'ACGT' ] ) ) >  1 ]
        with open('%s_positions_%s'%(outname_prefix, j), 'w' ) as out:
            out.write( '\n'.join( [ str( n+1 ) for n in nonshared_pos ] ) )
            if len(nonshared_pos)>0: 
                out.write('\n')
    
    out_q = multiprocessing.Queue()
    chunksize = int(math.ceil(reflen / float(nprocs)))
    procs = []
    
    for i in range(nprocs):
        p = multiprocessing.Process(
                target=subset_pos,
                args=(i, seqlist, masksites, i*chunksize, (i+1)*chunksize, out_q, outname_prefix))
        procs.append(p)
        p.start()
    
    
    # Wait for all worker processes to finish
    for p in procs:
        p.join()
    
    
    time_elapsed = datetime.datetime.now() - start
    time_elapsed = float(time_elapsed.seconds) + float(time_elapsed.microseconds/1e6)
    sys.stdout.write("Successfully wrote positions in %0.1f seconds.\n"%time_elapsed)
    sys.stdout.flush()
    start = datetime.datetime.now()
    
    cmd = 'cat %s_positions_* > %s_positions.txt'%(outname_prefix, outname_prefix)
    call(cmd, shell=True)
    cmd = 'rm %s_positions_*'%outname_prefix
    call(cmd, shell=True)
    
    # read back in non-shared positions
    posfile = '%s_positions.txt'%outname_prefix
    f = open(posfile, 'r')
    nonshared_pos = [int(l.strip())-1 for l in f]
    f.close()
    
    sys.stdout.write('Successfully obtained masked nonshared_diffs; there are %s of them.\n'%len( nonshared_pos ) )
    sys.stdout.flush()
    
    ## Sort out the SNPs
    for seq in seqlist:
        nonshared_bases = ''.join( seq.seq[ i ] for i in nonshared_pos )
        seq.seq._data = nonshared_bases
    SeqIO.write( seqlist , '%s_snps.fa'%outname_prefix, 'fasta' )
    sys.stderr.write('Successfully wrote snps fasta file.\n')
    sys.stdout.flush()
    
    ## Write the positions
    with open('%s_positions.txt'%outname_prefix, 'w' ) as out:
        out.write( '\n'.join( [ str( n+1 ) for n in nonshared_pos ] ) )
        out.write('\n')
    sys.stderr.write('Successfully wrote nonshared positions.\n')
    sys.stdout.flush()
    
    time_elapsed = datetime.datetime.now() - start
    time_elapsed = float(time_elapsed.seconds) + float(time_elapsed.microseconds/1e6)
    sys.stdout.write("Successfully completed alignment in %0.1f seconds.\n"%time_elapsed)
    sys.stdout.flush()