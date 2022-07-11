#!/usr/bin/python3
# -*- coding: utf-8 -*-
# distance.py
# D. Gibbon
# 2019-04-28

#==============================================================
# Input CSV file of property vectors for set to be displayed as dendrogram
# Property list format:

# NAME	
# 01-they-met-f:	1.27686245847	1.11473773411	1.00276754974	1.04728019047	...
# 02-I-will-be-:	2.60092490536	1.54899866894	1.53193261217	1.04048629977	...
# ...
#==============================================================
#==============================================================

# CONFIGURATION

def inputtextlines(filename):
	handle = open(filename,'r')
	linelist = handle.readlines()
	handle.close()
	return linelist

#==============================================================
# Graphviz network representation

def distnetworks(names, data, distances, netname, distmetric, mindistance, maxdistance):

	filename = netname + "-gv-" + distmetric

	# Select: dot, neato, fdp, twopi und circo

	d = Graph('D', filename=filename, engine='dot', format='png')
	d.attr('node', shape='ellipse', fontsize='12', size='6,6', rankdir='LR')

	dist_square = dist.squareform(distances)

	dist_list = dist_square.reshape(dist_square.shape[0] * dist_square.shape[1])

	dist_list = (dist_list - np.min(dist_list)) / (np.max(dist_list) - np.min(dist_list))

	ncount = len(dist_list)
	dist_square = dist_list.reshape(dist_square.shape)

	count = 0
	for i in range(0, len(names)-1):
		for j in range(i+1, len(names)):
			firstname = names[i]
			secondname = names[j]
			distance = dist_square[i][j]
			if distance >= mindistance and distance <= maxdistance:
				count += 1
				d.node(firstname)
				d.node(secondname)
				d.edge(firstname, secondname, label="%.3f"%distance)

	d.node(netname +"\n" + distmetric + ' distance metric\nn=%d/%d, %s min %s max'%(count,ncount,mindistance, maxdistance), shape='box')

	d.render(filename, view=False, format="png")

	return

#==============================================================

# Create dendrogram

def drawdendrogram(names, data, netname, dendroname, figparams, mindistance, maxdistance, showgraph):

	(figwidth, figheight, boxwidth, boxheight,
		halign, valign, orientation,
		netname, dendroname) = figparams

	data = np.array(data)
	# Length equalisation loop missing in earlier versions. Added on 2022-07-11, DG
	# Equalise lengths of data
	newsize = np.max( [ len(val) for val in data ] )
	valuelist = []
	for val in data:
		size = len(val)
		xloc = np.arange(size)
		new_xloc = np.linspace(0, size, newsize)
		new_data = np.interp(new_xloc, xloc, val)
		valuelist += [ new_data ]
	data = np.array(valuelist)

	#---------------------------------------------------------
	# Distance metrics
	"""
	# Exhaustive list for testing:
	distmetriclist = ['braycurtis', 'canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine',
		'dice', 'euclidean', 'hamming', 'jaccard', 'jensenshannon', 'kulsinski',
		'mahalanobis', 'matching', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean',
		'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
	"""

	# Short list: some of the other distance metrics yield zero or nan
	distmetriclist = ['canberra', 'chebyshev', 'cityblock', 'correlation', 'cosine', 'euclidean' ]

	print("Loop through similarity types: calculate distance table, create dendrogram")
	for distmetric in distmetriclist:

		print("----> Similarity type:", distmetric)
		distances = dist.pdist(data, metric=distmetric)

		#---------------------------------------------------------
		# GRAPHVIZ: similarity networks
		distnetworks(names, data, distances, netname, distmetric, mindistance, maxdistance)

		#---------------------------------------------------------
		methodlist = [ 'single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward' ]
		# Cycle through selected clustering methods
		for method in methodlist:

			print("===> Similarity type:", distmetric, "Clustering method:", method)
			dendrographname = "%-s-%s-%s"%(dendroname,distmetric,method)

			# Define figure
			fig = plt.figure(figsize=(figwidth, figheight))

			ax1 = fig.add_axes([halign, valign, boxwidth, boxheight])
			ax1.set_xlabel(dendrographname, fontsize=12)

			# Set plot font size for x-axis
#			plt.rc('xtick', labelsize=10)
			# rcParams.update({'font.size': 10})
			orientation = 'left'	# Change to 'right' or 'top' if leaf labels are cut off

			# Cluster calculation from distances and method
			Y1 = hy.linkage(distances, method=method)

			# Dendrogram creation from linkages and labels
			# cutoff = 0.3*np.max(Y1[:,2])
			hy.dendrogram(Y1,
				orientation=orientation,
				above_threshold_color='black', color_threshold=0,
				count_sort="False", distance_sort=False,
				labels=names, leaf_font_size=12)

			# 	plt.tight_layout()	# Seems incompatible with dendrogram figure

			plt.savefig(dendrographname + ".png")
			if showgraph:
				plt.show()

			plt.close(fig)	# Close each graph in loop after saving and displaying

	return

#==============================================================
#==============================================================
# Main caller

# Module import

try:
	import sys, re, time
	import numpy as np
	from matplotlib import rcParams
	import matplotlib.pyplot as plt
	import scipy.cluster.hierarchy as hy
	import scipy.spatial.distance as dist
	from graphviz import Graph
	import codecs
	encoding = "utf8"
	configfilename = "distance.conf"
	configfilename = re.sub(".py",".conf",sys.argv[0])
except:
	print("Import error."); sys.exit()

#-----------------------------------------------------

if True:

	if len(sys.argv) < 2:
		print("Error. Usage: %s <yournumvectorfile.csv>"%sys.argv[0])
		exit()

	csvfilename = sys.argv[1]
	filebase = re.sub(".csv", "", re.sub(".*/", "", csvfilename))
	dendroname = "DENDRO/"+filebase + "-dendro"
	netname = "GRAPHVIZ/" + filebase + "-network"

	from numdistnetdendro_conf import *
	"""
	#---------------------------------------------------------------
	# Assign configuration parameters to variables (see .conf file)
	conflines = [
		re.sub(" ","",line.rstrip()).split("=")
		for line in inputtextlines(configfilename)
		if line[0] != '#'
		]
	confdict = { line[0] : line[1] for line in conflines if line != [] and line[0] != '' }

	showgraph = False if confdict['showgraph'] == 'False' else True
	figwidth = float(confdict['figwidth'])
	figheight = float(confdict['figheight'])
	boxwidth = float(confdict['boxwidth'])
	boxheight = float(confdict['boxheight'])
	halign = float(confdict['halign'])
	valign = float(confdict['valign'])
	orientation = confdict['orientation']
	mindistance = float(confdict['mindistance'])
	maxdistance = float(confdict['maxdistance'])
	"""
	# Dendrogram parameters
	figparams = (figwidth, figheight, boxwidth, boxheight,
		halign, valign, orientation, netname, dendroname)

if False:
	print("Configuration error."); exit()

# INPUT CSV TABLE AS NAME LIST AND VALUELIST, GENERATE DENDROGRAM

if True:
	# Preprocess CSV strings
	csvlist = inputtextlines(csvfilename)
	csvlist = [re.sub(" ","",line.rstrip()).split(",") for line in csvlist ]

	# Extract name list and vector list from propvectors
	names = [ row[0] for row in csvlist ]
	data = [ [ float(chars) for chars in row[1:]] for row in csvlist ]

	drawdendrogram(names, data, netname, dendroname, figparams, mindistance, maxdistance, showgraph)

if False:
	print("Dendrogram creation error."); sys.exit()


#EOF===========================================================
