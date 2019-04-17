#!/usr/bin/env python

#############################################################################
# clusterCreator.py
# takes a list of pairwise differences and clusters related samples
#############################################################################


from optparse import OptionParser
import os, sys, subprocess
import networkx as nx

usage = '''
Dependencies - networkx: pip install networkx --user

clusterCreator instructions
Usage: clusterCreator.py [options] node_list.txt edge_list.txt output_file

node_list.txt should be a list of nodes, 1 on each line, no header
edge_list.txt should be a list of edges, in tab separated format, no header: node1 node2 snvs

Example:

./clusterCreator -s 10 -m 2 -r C00000001,C00000014 nodes.txt edges.txt out.txt

Options
-s, --snvs			snv cutoff value for clustering, must be within this limit to be placed in a cluster [default = 100]
-m, --mincc			minimum sized connected component to keep, use if only want to keep clusters above a certain size [default = 1]
-r, --remove		comma separated list of nodes to remove, no spaces after comma
-f, --removefile	file containing list of nodes to remove
'''


if __name__=="__main__":
	

	### PROCESS INPUT OPTIONS AND ARGUMENTS ###	
	parser = OptionParser()
	parser.add_option( "-s", "--snvs", action = "store", type="int", dest = "snv_threshold", default = 100 )
	parser.add_option( "-m", "--mincc", action = "store", type="int", dest = "mincc", default = 1 )
	parser.add_option( "-r", "--remove", action = "store", dest = "remove_nodes", default = "" )
	parser.add_option( "-f", "--removefile", action = "store", dest = "remove_file", default = "" )
	opts, args = parser.parse_args()
	
	#check have correct number of arguments
	if len(args) == 3:
		node_file, edge_file, outfile = args
		snv_threshold = opts.snv_threshold
		mincc = opts.mincc
		remove_nodes = opts.remove_nodes.split(',')
		remove_file = opts.remove_file
		#read in remove_file
		if remove_file != "":
			remove_nodes = []
			f = open(remove_file, 'r')
			for l in f:
				if not l.startswith('#'):
					remove_nodes.append(l.strip())
	else:
		sys.stdout.write(usage)
		sys.exit(1)


	## echo out removed nodes
	if len(remove_nodes)>0:
		sys.stdout.write("Removing excluded nodes\n")
		sys.stdout.write("\n".join(remove_nodes))
		sys.stdout.flush()
	
	### READ IN FILES ###	
	#read in nodes
	f = open(node_file, 'r')
	nodes = [l.strip() for l in f]
	f.close()

	#read in edges
	f = open(edge_file, 'r')
	edges = [l.strip().split()[0:2] for l in f if int(l.strip().split()[2])<=snv_threshold]
	f.close()

	### CREATE GRAPH ###	
	#create a networkx graph and add the nodes and edges
	G = nx.Graph()
	G.add_nodes_from(nodes)
	G.add_edges_from(edges)
	G.remove_nodes_from(remove_nodes)

	#find the connected components, i.e. the clusters
	CC = nx.connected_components(G)

	#open output file
	o = open(outfile, 'w')

	#add header, where cluster_number is an identifier for the cluster
	o.write('cluster_number\tid\n')

	#i is the cluster_number
	i = 0
	for cc in CC:
		if len(cc)>=mincc:
			i +=1
			for n in cc:
				o.write('%s\t%s\n'%(i, n))

	o.close()
