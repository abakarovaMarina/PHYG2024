# -*- coding: Utf-8 -*-

# developped by Aubin Fleiss
# contact : aubin.fleiss@gmail.com

from random import randint, choice, randrange, shuffle, seed, random as rd
import numpy as np
import copy
from ete2 import *
from itertools import islice

# Note : fusion can happen only if the number of chromosoms >=2 i have to implement this
# Note : same for transocation

def poisson(poisson_mean):
	"""This function return a value according to a Poisson distribution
	centered around the Poisson_mean value"""
	val = 0
	while val == 0:
		val = int(np.random.poisson(poisson_mean, 1))
	return val


def choose_coordinates(genome):
	"""this function chooses a random position in the genome"""

	# choose chromosome : adjust probability based on chromosome length
	chroms=genome.keys()
	nbGenes=sum([len(genome[chrom]) for chrom in genome])
	p=[float(len(genome[chrom]))/float(nbGenes) for chrom in sorted(genome.keys())]
	chromChoice = int(np.random.choice(sorted(genome.keys()), 1, p=p))

	# choose gene in previously selected chromosome
	geneChoice=choice(range(len(genome[chromChoice])))

	return((chromChoice,geneChoice))


def reverse(l):
	return [ -1*val for val in l[::-1] ]

def reciprocal_translocation(genome):
	#xprint("reciprocal_translocation")

	# TEMPORARY FIX
	if len(genome.keys()) < 2:
		return(genome)

	#choose two positions
	while True:
		coord1 = choose_coordinates(genome)
		cen1   = genome[coord1[0]].index(0)
		if cen1!=coord1[1] and coord1[1]!=0 and coord1[1]!=len(genome[coord1[0]])-1:
			break

	coord2 = coord1
	while True:
		coord2 = choose_coordinates(genome)
		if coord1[0]==coord2[0]: # pick second position from a different chromosome
			continue
		cen2 = genome[coord2[0]].index(0)
		if cen2!=coord2[1] and coord2[1]!=0 and coord2[1]!=len(genome[coord2[0]])-1:
			break

	#This is where the translocation is implemented
	len1= len(genome[coord1[0]])
	len2= len(genome[coord2[0]])

	cut1 = coord1[1]
	cut2 = coord2[1]

	if cut1<cen1 and cut2<cen2:
		##xprint("case1")
		ch1 = genome[coord2[0]][0:cut2] + genome[coord1[0]][cut1:]
		ch2 = genome[coord1[0]][0:cut1] + genome[coord2[0]][cut2:]
	elif cut1<cen1 and cut2>cen2:
		##xprint("case2")
		cutted1 = genome[coord1[0]][0:cut1]
		cutted2 = genome[coord2[0]][cut2+1:]
		ch1 = reverse(cutted2) + genome[coord1[0]][cut1:]
		ch2 = genome[coord2[0]][0:cut2+1] + reverse(cutted1)
	elif cut1>cen1 and cut2<cen2:
		##xprint("case3")
		cutted1 = genome[coord1[0]][cut1+1:]
		cutted2 = genome[coord2[0]][0:cut2]
		ch1 = genome[coord1[0]][0:cut1+1] + reverse(cutted2)
		ch2 = reverse(cutted1) + genome[coord2[0]][cut2:]
	elif cut1>cen1 and cut2>cen2:
		##xprint("case4")
		ch1 = genome[coord1[0]][0:cut1+1] + genome[coord2[0]][cut2+1:]
		ch2 = genome[coord2[0]][0:cut2+1] + genome[coord1[0]][cut1+1:]

	genome[coord1[0]] = ch1
	genome[coord2[0]] = ch2

	return(genome)


def aux1(genome, coord, interval):
	values=[]
	nb=0
	if len(genome[coord[0]])<interval:
		return False
	while len(values)<interval:
		try:
			if coord[1]+nb<0:
				raise
			else:
				genome[coord[0]][coord[1]+nb]
				values.append(nb)
		except:
			pass

		if nb<=0:
			nb=-1*nb+1
		else:
			nb=-1*nb

		values=sorted(values)
	return(sorted([coord[1]+max(values),coord[1]+min(values)]))


def duplication(genome,mean_dup_len,allow_mirror_dup=False):
	#xprint("duplication")

	while True:
		coord=choose_coordinates(genome)
		interval = poisson(mean_dup_len)
		ex = aux1(genome, coord, interval)
		if ex!=False:
			break

	temp  = [ genome[coord[0]][i] for i in range(min(ex),max(ex)+1)]
	temp1 = list(temp)
	temp2 = list(temp)

	remove_duplicated_centromere=False
	if 0 in temp:
		# if we try to duplicate the centromere, keep in mind to remove one of them
		remove_duplicated_centromere=True

	if allow_mirror_dup:
		mirror = choice([True,False])
		if mirror:
			temp2 = [-1*elt for elt in temp2]
			temp2.reverse()

	# this is where the duplication is implemented
	genome[coord[0]]=genome[coord[0]][0:ex[0]]+temp1+temp2+genome[coord[0]][ex[1]+1:]

	if remove_duplicated_centromere:
		centromere_indexes = [i for i,x in enumerate(genome[coord[0]]) if x == 0]
		removed_centromere = choice(centromere_indexes)
		genome[coord[0]] = [gene for i,gene in enumerate(genome[coord[0]]) if i!=removed_centromere ]

	return(genome)


def fusion(genome):
	#xprint("fusion")

	# TEMPORARY FIX
	if len(genome.keys()) < 2:
		return(genome)

	ch1 = choice(genome.keys())
	ch2 = ch1
	while ch2 == ch1:
		ch2=choice(genome.keys())

	# delete one of the two centromeres here
	acentromeric_chrom = choice((ch1,ch2))
	genome[acentromeric_chrom] = [gene for gene in genome[acentromeric_chrom] if gene!=0]

	ext1=choice(("begin","end"))
	ext2=choice(("begin","end"))

	if ext1=="end":
		if ext2=="begin":
			##xprint("case1")
			new_chrom = genome[ch1]+genome[ch2]
		if ext2=="end":
			##xprint("case2")
			new_chrom = genome[ch1]+reverse(genome[ch2])

	if ext1=="begin":
		if ext2=="begin":
			##xprint("case3")
			new_chrom= reverse(genome[ch2])+genome[ch1]
		if ext2=="end":
			##xprint('case4')
			new_chrom= genome[ch2]+genome[ch1]

	del genome[ch1]
	del genome[ch2]
	genome[min(ch1,ch2)] = new_chrom

	return(genome)


def wgd(genome,averageDeletionRateWGD):
	#xprint("whole_genome_duplication")
	genome2=dict()
	for k in genome.keys():
		genome2[k]=[elt for elt in genome[k]]

	for i in sorted(genome.keys()):
		for j in range(len(genome[i])):
			keep_both = np.random.choice((True,False), 1, p=(1-averageDeletionRateWGD,averageDeletionRateWGD))
			if genome[i][j]==0:
				keep_both=True

			if keep_both:
				pass
			else:
				version = choice((1,2))
				if version==1:
					genome[i][j]  = None
				else:
					genome2[i][j] = None

	for i in sorted(genome.keys()):
		temp = []
		for gene in genome[i]:
			if gene is not None:
				temp.append(int(gene))
		genome[i]=temp

	for i in sorted(genome2.keys()):
		temp = []
		for gene in genome2[i]:
			if gene is not None:
				temp.append(int(gene))
		genome2[i]=temp

	# NB
	# centromeres are not renamed, whatever happens
	# genes are renamed when they form paralogs
	# otherwise they keep their name even if they are
	# not located on their original chromosome

	assembly = [genome[ch] for ch in sorted(genome.keys())]+ [genome2[ch] for ch in sorted(genome2.keys())]
	i=1
	genome=dict()

	for ch in assembly :
		genome[i]=ch
		i+=1


	return(genome)

def print_genome(genome):
	for chrom in genome.keys():
		print(chrom),
		print(genome[chrom])

# ---------------- tree generator -------------------

def make_tree(nbSpecies):

	# Generate tree Raphael style
	t=Tree()
	for i in range(nbSpecies-1):
		node=t.get_tree_root()
		while not node.is_leaf():
			node = choice(node.get_children())
		node.add_child()
		node.add_child()
	for node in t.traverse():
		node.dist=rd()

	# generate tree ete2 style
	t=Tree()
	t.populate(nbSpecies, random_branches=True)
	t.ladderize()

	# rename species
	spName=0
	anName=0
	for node in t.traverse("preorder"):
		if node.is_leaf():
			spName+=1
			node.name = "s"+"0"*(4-len(str(spName))-1)+str(spName)
		elif not node.is_leaf():
			anName+=1
			node.name= "a"+"0"*(4-len(str(anName))-1)+str(anName)

	# expand leaves to the same level
	leaves = [node for node in t.traverse() if node.is_leaf()]
	max_dist=max([leaf.get_distance(t.get_tree_root().name) for leaf in leaves])

	for leaf in leaves:
		parent = leaf.up
		leaf.dist = max_dist-parent.get_distance(t.get_tree_root().name)

	for node in t.traverse():
		node.dist= node.dist/max_dist

	return(t)



def random_chunk(li, min_chunk=1, max_chunk=3):
	it = iter(li)
	while True:
		nxt = list(islice(it,randint(min_chunk,max_chunk)))
		if nxt:
			yield nxt
		else:
			break

def random_split(seq, n):
	rnd_bools = random.sample((0,) * n + (1,) * (len(seq) - n), len(seq))
	left_right = ([], [])
	for b, x in zip(rnd_bools, seq):
		left_right[b].append(x)
	return left_right

def make_root_genome(nbChroms,nbGenes):

	if nbGenes/nbChroms<10:
		#xprint("The genome must have a mean of at least 10 genes per chromosome")
		exit()
	while True:
		temp=list(range(1,nbGenes+1))
		genes_per_chrom = list(np.random.normal(nbGenes/nbChroms,0.4*nbGenes/nbChroms,nbChroms-1))
		genes_per_chrom = [ int(elt) for elt in genes_per_chrom ]
		genes_per_chrom.append(nbGenes-sum(genes_per_chrom))
		if False not in [elt>5 for elt in genes_per_chrom]:
			break

	genome={i+1:range(sum(genes_per_chrom[0:i])+1,sum(genes_per_chrom[0:i])+genes_per_chrom[i]+1) for i in range(len(genes_per_chrom))}

	# add centromeres
	for ch in sorted(genome.keys()):
		while True:
			centro_pos = choice(range(len(genome[ch])))
			if centro_pos!=0 and genome[ch][centro_pos]!=genome[ch][-1]:
				break
		temp=list(genome[ch])
		temp.insert(centro_pos,0)
		genome[ch]=temp

	return(genome)
