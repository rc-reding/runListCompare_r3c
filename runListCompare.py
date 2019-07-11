#!/usr/bin/env python

## takes a single arguement which is the path to an ini file

from subprocess import call
from subprocess import check_output
import sys, os, gzip
from configparser import ConfigParser
from Bio import SeqIO

import treeswift


### import options from ini file ###
ini_file = sys.argv[1]

parser = ConfigParser()
parser.read(ini_file)

#import reference section
refpath = parser.get('ref', 'refpath')
reflen = int(parser.get('ref', 'reflen'))
ref = parser.get('ref', 'ref')
refname = parser.get('ref', 'refname')
maskfile = parser.get('ref', 'maskfile')

#import samples section
seqlist = parser.get('samples', 'seqlist')

#import options
perACGT_cutoff = float(parser.get('options', 'perACGT_cutoff'))
nprocs = int(parser.get('options', 'nprocs'))
cluster_snp = int(parser.get('options', 'cluster_snp'))
varsite_keep = float(parser.get('options', 'varsite_keep'))
seq_keep = float(parser.get('options', 'seq_keep'))
align_n = int(parser.get('options', 'align_n'))

#output options
output_stem = parser.get('output', 'output_stem')
round_dp = int(parser.get('output', 'round_dp'))
draw_cf = parser.getboolean('output', 'draw_cf')
use_pypy = parser.getboolean('output', 'use_pypy')
use_iqtree = parser.getboolean('output', 'use_iqtree')


### end inputs ###

sys.stdout.write('\nChecking percentage ACGT in samples\n')

def getPct(fapath, reflen):
    if os.path.exists( fapath):
        sys.stdout.write('Checking percentage ACGT in: %s\n'%fapath)
        sys.stdout.flush()
        with gzip.open(fapath, 'rt') as fa_fh:
            record = list(SeqIO.parse(fa_fh, 'fasta'))[0]
        # May want to consider possiblity of multiple contigs here
        # Deleted some code which did?
        len_acgt = len([b for b in record.seq if b in 'ACGT'])
        return len_acgt/reflen
    else:
        sys.stderr.write(fapath+' does not exist, skipping...')
        return None


## check if directories exist, and if not make
# output stem directory
if not os.path.isdir(output_stem):
    os.mkdir(output_stem)
    
# copy ini file to output directory
cmd = 'cp %s %s'%(ini_file, output_stem)
call(cmd.split())

# reject directory
if not os.path.isdir('%s/reject'%output_stem):
    os.mkdir('%s/reject'%output_stem)

# cluster directory
if not os.path.isdir('%s/cluster'%output_stem):
    os.mkdir('%s/cluster'%output_stem)

# cluster_ml directory
if not os.path.isdir('%s/cluster_ml'%output_stem):
    os.mkdir('%s/cluster_ml'%output_stem)

## read in list of fasta files
f = open(seqlist, 'r')
seqlist = [l.strip().split() for l in f if l.strip()]

#check all sequences are of sufficient length and write pctACGT
f = open('%s/pctACGT.txt'%output_stem, 'w')
clean_seqlist = []
for s in seqlist:
    l = getPct(s[1], reflen)
    if l:
        f.write('%s\t%s\t%0.3f\n'%(s[0], s[1], l))
        if l>perACGT_cutoff:
            clean_seqlist.append(s)

f.close()


#write clean seqlist
temp_seqlist = '%s/clean_seqlist.txt'%output_stem
f = open(temp_seqlist, 'w')
for s in clean_seqlist:
    f.write('%s\t%s\n'%(s[0], s[1]))
f.close()

## run initiall all vs all comparison
#generate alignment
pypy_stem = ''
if not use_pypy:
    pypy_stem = 'python '


sys.stdout.write('\nGenerate all vs all alignment using %s cores\n'%nprocs)
if maskfile:
    sys.stdout.write('Proceeding with mask file\n')
    cmd = '%smtAlign.py -p %s -m %s %s %s %s/align'%(pypy_stem, nprocs, maskfile, temp_seqlist, refpath, output_stem)
else:
    sys.stdout.write('Proceeding without mask file\n')
    cmd = '%smtAlign.py -p %s %s %s %s/align'%(pypy_stem, nprocs, temp_seqlist, refpath, output_stem)

sys.stdout.write(cmd+'\n')
sys.stdout.flush()
call(cmd.split())

sys.stdout.write('\nGenerate all vs all distance matrix %s cores\n'%nprocs)
cmd = 'python getDist.py -p %s %s/align_snps.fa %s/align-compare'%(nprocs, output_stem, output_stem)
sys.stdout.write(cmd+'\n')
sys.stdout.flush()
call(cmd.split())


#merge alignment files
cmd = 'cat %s/align-compare_* >%s/align-compare.txt'%(output_stem, output_stem)
call(cmd, shell=True)

cmd = 'rm %s/align-compare_*'%output_stem
call(cmd, shell=True)


## cluster into groups within 100 SNVs
# generate node file
cmd = 'cut -f 1 %s/clean_seqlist.txt > %s/initial_nodes.txt'%(output_stem, output_stem)
call(cmd, shell=True)

sys.stdout.write('\nGenerate clusters and clean alignments\n')
## clear reject folder
cmd = 'rm -f %s/reject/*'%output_stem
call(cmd, shell=True)   

def run_cluster(maskfile, cluster_snp, exclude, output_stem, nprocs, refpath, seq_keep, varsite_keep):
    # generate cluster file 
    if exclude == '':
        r = ''
    else:
        r = '-r %s'%exclude
    
    cmd = 'python clusterCreator.py -s %s %s %s/initial_nodes.txt %s/align-compare.txt %s/clusters.txt'%(cluster_snp, r, output_stem, output_stem, output_stem)
    sys.stdout.write(cmd+'\n')
    sys.stdout.flush()
    call(cmd.split())
    
    ## clean cluster directories
    cmd = 'rm -f %s/cluster/*'%output_stem
    call(cmd, shell=True)
    
    cmd = 'rm -f %s/cluster_ml/*'%output_stem
    call(cmd, shell=True)
    
    
    ## run cluster alignment and cleaning
    if maskfile:
        cmd = 'python getClusterAlign.py -p %s -s %s -v %s -n %s -m %s %s/clean_seqlist.txt %s/clusters.txt %s %s'%(nprocs, seq_keep, varsite_keep, align_n, maskfile, output_stem, output_stem, refpath, output_stem)
    else:
        cmd = 'python getClusterAlign.py -p %s -s %s -v %s -n %s %s/clean_seqlist.txt %s/clusters.txt %s %s'%(nprocs, seq_keep, varsite_keep, align_n, output_stem, output_stem, refpath, output_stem)
    sys.stdout.write(cmd+'\n')
    sys.stdout.flush()
    call(cmd.split())






### need to search for files that cause rejection and re-run 
# set up first run
run_cluster(maskfile, cluster_snp, '', output_stem, nprocs, refpath, seq_keep, varsite_keep)
# remove rejected nodes file
rejectedNodesFile = '%s/rejected_nodes.txt'%output_stem
if os.path.exists(rejectedNodesFile):
    os.remove(rejectedNodesFile)


# check for any reject files as loop
iteration = 0
while True:
    iteration +=1
    cmd = 'ls -l %s/cluster/*reject* 2> /dev/null | wc -l'%output_stem
    r = check_output(cmd, shell=True, text=True)
    r = int(r.strip())
    if r == 0:
        break
    
    # get list of rejected nodes
    cmd = 'cat %s/cluster/*reject*'%output_stem
    exclude = check_output(cmd, shell=True, text=True)
    # save these to file too
    cmd = 'cat %s/cluster/*reject* >> %s/rejected_nodes.txt 2>/dev/null'%(output_stem, output_stem)
    call(cmd, shell=True)
    
    #get the numbers of any rejected clusters - cannot use cluster_ in filename - copy and rename
    cmd = 'ls %s/cluster/*reject* 2> /dev/null'%output_stem
    reject = check_output(cmd, shell=True, text=True).split()
    for r in reject:
        r = r.split('/cluster/cluster_')[1].split('_')[0]
        cmd = 'cp %s/cluster/*_%s_* %s/reject/'%(output_stem, r, output_stem)
        call(cmd, shell=True)
    
    cmd = "rename 's/(.*)/$1_run_%s/' %s/reject/*"%(iteration, output_stem)
    call(cmd, shell=True)
    
    #read in previous exclusions and add these to string
    e = open('%s/rejected_nodes.txt'%output_stem, 'r')
    exclude_old = ','.join(set([l.strip().split()[0] for l in e]))
    e.close()
    
    exclude = ','.join([e.split()[0] for e in exclude.strip().split('\n')])
    
    #read in old excludes
    exclude = ",".join([exclude_old,exclude])
    
    # loop around
    sys.stdout.write('Repeating analysis with Reject: %s\n'%exclude)
    sys.stdout.flush()
    run_cluster(maskfile, cluster_snp, exclude, output_stem, nprocs, refpath, seq_keep, varsite_keep)

### run ML comparison

sys.stdout.write('\nGenerate ML trees\n')
sys.stdout.flush()

tree = ''
if draw_cf: tree = '-c '
iqtree = ''
if use_iqtree: iqtree = '-q '
cmd = 'python getClusterML.py %s%s-p %s -r %s %s %s'%(tree, iqtree, nprocs, round_dp, refpath, output_stem)
sys.stdout.write(cmd+'\n')
sys.stdout.flush()
call(cmd, shell=True)

# write ML distances file


cmd = 'ls %s/cluster_ml/*scale*'%output_stem  # Too low a SNP threshold and this fails
files = check_output(cmd, shell=True, text=True)
files = files.split()

outfile_ml = '%s/ML_distances.txt'%output_stem
w_ml = open(outfile_ml, 'w')

if draw_cf:
    outfile_cf = '%s/CF_distances.txt'%output_stem
    w_cf = open(outfile_cf, 'w')

for f in files:
    tr = treeswift.read_tree_newick(f)
    dist_mat = tr.distance_matrix(leaf_labels=True)
    dists = {(k, kk): vv  # Make flat tuple-keyed dict from nested dict of pairwise dists
             for k, v in list(dist_mat.items())
             for (kk, vv) in list(v.items())}
    for pair in list(dists.keys()):
        if 'phyml_tree_scaled.tree' in f:
            w_ml.write('%s\t%s\t%s\n'%(pair[0], pair[1], dists[pair]))
        if 'cf_scaled.tree' in f and draw_cf:
            w_cf.write('%s\t%s\t%s\n'%(pair[0], pair[1], dists[pair]))

w_ml.close()
if draw_cf: w_cf.close()

sys.stdout.write('Complete\n\n')
sys.stdout.flush()

## final outputs as follows in cd-link
# ML_distances.txt - ML distances
# CF_distances.txt - CF distances
# align-compare.txt - raw pairwise distances, includes excluded nodes
# rejected_nodes.txt - list of excluded nodes based on local alignment cleaning
# clusters.txt - final list of clusters
# cluster folder - all alignments prior to and post cleaning, but excluding nodes
# cluster_ml - ml files for each cluster

