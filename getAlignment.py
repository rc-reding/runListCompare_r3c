#!/usr/bin/env python

## Script generates an alignment of variable sites, a list of positions
## Allows for masking of given elements of genome

## uses pypy rather than python 

import sys, re, os, gzip, time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

from optparse import OptionParser

usage = '''
getAlignment instructions
Usage: getAlignment.py [options] list_of_fasta.txt ref.fa output_prefix

list_of_fasta.txt: tab delimited 1) name, 2) fasta path and filename
ref.fa: reference in fasta format
output_prefix: output folder

Options
-m, --mask  optional mask file for mobile elements, 
            tab delimited - start and end coordinates, counting from one, not zero

Example:

./getAlignment.py -m mask.txt list_of_fasta.txt ref.fa output_prefix


'''


if __name__=='__main__':

    start = time.perf_counter()

    ### PROCESS INPUT OPTIONS AND ARGUMENTS ### 
    parser = OptionParser()
    parser.add_option( '-m', '--mask', action = 'store', type='string', dest = 'maskfile', default = '' )
    opts, args = parser.parse_args()
    
    #check have correct number of arguments
    if len(args) == 3:
        listoffasta, refpath, outname_prefix = args
        maskfile = opts.maskfile
    else:
        sys.stdout.write(usage)
        sys.exit(1)
    
    #read in reference
    #refseq = SeqIO.read( refpath, 'fasta' ).seq
    # refseq = "".join([s.seq._data for s in SeqIO.parse( refpath, 'fasta' )])
    refseq = next(SeqIO.parse( refpath, 'fasta'))
    reflen = len(refseq.seq)

    sys.stderr.write('Successfully read in reference.\n')    
    
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
    with open(listoffasta) as fp:
        for line in fp:
            nicename, fapath = line.strip().split()
            if (os.path.exists( fapath)):
                sys.stdout.write('%s\n'%fapath)
                with gzip.open(fapath, 'rt') as fa_fh:
                    fa = list(SeqIO.parse(fa_fh, 'fasta'))[0]
                    fa.id = nicename
                    fa.description = ''
                seqlist.append(fa)
            else:
                sys.stderr.write(fapath+' does not exist, skipping...')
    
    time_elapsed = time.perf_counter() - start
    sys.stdout.write("Successfully read in sequences in %s seconds.\n"%time_elapsed)
    start = time.perf_counter()
    
    ## find nonshared positions after masking
    seq_generator = zip( masksites, *seqlist  )
    nonshared_pos =[ i for ( i, a ) in enumerate( seq_generator ) 
                   if a[0]==1 and len( set( [ ai for ai in a[1:] if ai in 'ACGT' ] ) ) >  1 ]
    sys.stdout.write('Successfully obtained masked nonshared_diffs; there are %s of them.\n'%len( nonshared_pos ) )
    
    ## Sort out the SNPs
    for seq in seqlist:
        nonshared_bases = ''.join( seq.seq[ i ] for i in nonshared_pos )
        seq.seq._data = nonshared_bases
    SeqIO.write( seqlist , '%s_snps.fa'%outname_prefix, 'fasta' )
    sys.stderr.write('Successfully wrote snps fasta file.\n')
    
    ## Write the positions
    with open('%s_positions.txt'%outname_prefix, 'w' ) as out:
        out.write( '\n'.join( [ str( n+1 ) for n in nonshared_pos ] ) )
        out.write('\n')
    sys.stderr.write('Successfully wrote nonshared positions.\n')
    
    time_elapsed = time.perf_counter() - start
    sys.stdout.write("Successfully completed alignment in %s seconds.\n"%time_elapsed)

