# -*- coding: Utf-8 -*-

# developped by Aubin Fleiss
# contact : aubin.fleiss@gmail.com

import random as rd
import numpy as np
import copy
import sys

""" This function return a non-zero value according to a Poisson distribution
with the lambda parameter provided in input """
def poisson(lam):
    val = 0
    while val == 0:
        val = int(np.random.poisson(lam, 1))
    return val

""" This function takes in input a genome and returns a pair where
the first value is a chromosome identifier and the second value is a random
index within the list of genes of the choosen chromosome """
def choose_coordinates(genome):
    # choose chromosome : adjust probability based on chromosome length
    chroms=genome.keys()
    nbGenes=sum([len(genome[chrom]) for chrom in genome])
    p=[float(len(genome[chrom]))/float(nbGenes) for chrom in sorted(genome.keys())]
    chromChoice = int(np.random.choice(sorted(genome.keys()), 1, p=p))
    # choose gene in previously selected chromosome
    geneChoice = rd.choice(range(len(genome[chromChoice])))
    return((chromChoice,geneChoice))

""" This function takes in input a genome and prints it to the standard output """
def print_genome(genome):
    for chrom in genome.keys():
        print(chrom),
        print(genome[chrom])


def inversion(genome,mean_inv_len):
    """ EXERCISE 2 - Inversion function """
    return(genome)


def deletion(genome,mean_del_len):
    """ EXERCISE 2 - Deletion function """
    return(genome)


def fission(genome):
    """ EXERCISE 2 - Fission function """
    return(genome)


def main( argv=None ):

    # definition of a test genome
    genome = {}
    genome[1]=chr 2range(1,22)
    genome[1].insert(5,0)
    genome[2]=range(22,54)
    genome[2].insert(10,0)
    genome[3]=range(54,76)
    genome[3].insert(7,0)
    genome[4]=range(76,101)
    genome[4].insert(11,0)

    print("### ancestor genome ###")
    print_genome(genome)

    print("### genome after an inversion ###")
    genome = inversion(genome,5)
    print_genome(genome)

    print("### genome after a deletion ###")
    genome=deletion(genome,1)
    print_genome(genome)

    print("### genome fission ###")
    genome=fission(genome)
    print_genome(genome)

    return 0


if __name__=="__main__":
    sys.exit(main())
