# module_dendrogram.py
# D. Gibbon
# Created 2021-08-16
# Modified 2022-04-18
# Formant dendrogram for rfa.py

#============================================================

import matplotlib.pyplot as plt
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hy
import warnings
warnings.filterwarnings("ignore")

#============================================================
# Draw dendrogram of distances between peak frequencies, i.e. a list of length 1 vectors :)

def drawdendrogram(pltobj1, positionlist, boxwidth, boxheight, halign, valign, fontsize):

	poslist = list(positionlist)		# convert zip object to list of pairs
	namelist, valuelist  = zip(*poslist)
	namelistlen = len(namelist)		# get length of list of names

	namestringlist = [
		"f:%.3f\nT:%.3f"%(n, v)
		for n, v in zip(namelist, valuelist) ]
	valuelist = [ [x] for x in namelist]	# 2D array]

	# Distance metric
	metric = "cityblock"
	distances = dist.pdist(valuelist, metric=metric)

	# Distance dimension reduction (clustering)
	method = "ward"
	Y1 = hy.linkage(distances, method=method, optimal_ordering=True)

	# Dendrogram figure sizes
	dendrofontsize = 6
	orientation = "top"

	# Dendrogram figure - new overlaid box
	pltobj2 = plt.axes([halign,valign,boxwidth,boxheight])
	pltobj2.set_xticks([])
	pltobj2.set_yticks([])
	pltobj2.spines["top"].set_visible(False)
	pltobj2.spines["bottom"].set_visible(False)
	pltobj2.spines["left"].set_visible(False)
	pltobj2.spines["right"].set_visible(False)

	#	cutoff = 0.3*np.max(Y1[:,2])
	#	hy.leaves_list(Y1)	# To be used later

	# Dendrogram creation
	hy.dendrogram(
		Y1, ax=pltobj2, orientation=orientation,
		above_threshold_color='black',
		color_threshold=0,
		count_sort="False", distance_sort=False,
		labels=namestringlist, leaf_font_size=dendrofontsize)

	return

# EOF
