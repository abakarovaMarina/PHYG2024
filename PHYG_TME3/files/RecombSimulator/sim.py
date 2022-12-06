# -*- coding: Utf-8 -*-

# developped by Aubin Fleiss
# contact : aubin.fleiss@gmail.com

# import standard modules
from multiprocessing import Lock, Process, Queue, current_process, Pool, cpu_count
from random import randint, choice, randrange, shuffle
from os import getpid, makedirs, getcwd
from os.path import exists
import time
import numpy as np
from shutil import rmtree
import datetime
import time
import argparse

from os import listdir
from os.path import isfile, join

# import other modules
from ete2 import *

### todo list :
# - add RAM saturation warning for big simulations
# - generate stats output in pdf
# - implement log mode (where every rearrangement is written)

#======================================== simulation parameters ========================================

nbSimulations = 3 # number of simulations to perform
nbSpecies     = 10   # number of species in simulated tree
nbChroms      = 8   # number of chromosomes in ancestor
nbGenes       = 5000 # number of genes in ancestor
E = 500 # average number of events from the root to any species (~500 for yeasts, 1000 for vertebrates)
threadNum = int(cpu_count())

parameters={}
parameters["mean_dup_len"]     = 5 # mean duplication length (in genes)
parameters["mean_del_len"]     = 1 # mean deletion length (in genes)
parameters["mean_inv_len"]     = 5 # mean inversion length (in genes)
parameters["allow_mirror_dup"] = False # allow mirror duplications to happen, not only tandem duplications
parameters["averageDeletionRateWGD"] = 0.8 # probability to delete one of two copies of a duplicated gene after a WGD

# Probability of each kind of event (it is OK if they do not sum up to 1, they will be divided by their sum)
probas                  = {}
probas["Inversion"]     = 0.6     # probability of an inversion to happen
probas["Translocation"] = 0.2979  # probability of a reciprocal translocation to happen
probas["Duplication"]   = 0.05    # probability of a duplication to happen
probas["Deletion"]      = 0.05    # probability of a deletion to happen
probas["Fusion"]        = 0.001   # probability of a fusion to happen
probas["Fission"]       = 0.001   # probability of a fission to happen
probas["WGD"]           = 0.0001  # probability of a whole genome duplication to happen

minimalEventsOnBranch = 1 # minimal number of events on a branch # best not to modify

save_directory="results" # where the results will be saved => do not modify if possible
Chronicle_format=True # generate Chronicle files or not ? => do not modify

#======================================== get custom user parameters ======================================

parser = argparse.ArgumentParser()
parser.add_argument('--genes',   dest='geneNum',    type=int, required=True, help='number of genes')
parser.add_argument('--chr',     dest='chrNum',     type=int, required=True, help='number of chromosomes')
parser.add_argument('--species', dest='speciesNum', type=int, required=True, help='number of species')
parser.add_argument('--events',  dest='eventNum',   type=int, required=True, help='number of events from root to species')
parser.add_argument('--nsim',    dest='simNum',     type=int, required=True, help='number of simulations')
parser.add_argument('--cpu',     dest='threadNum',  type=int, help='number of CPUs to use')
args = parser.parse_args()

nbSimulations = args.simNum     # number of simulations to perform
nbSpecies     = args.speciesNum # number of species in simulated tree
nbChroms      = args.chrNum     # number of chromosomes in ancestor
nbGenes       = args.geneNum    # number of genes in ancestor
E = args.eventNum

if args.threadNum is not None:
    threadNum = min( threadNum, int(args.threadNum) )
if threadNum < 1:
    threadNum = 1

#======================================== simulation functions ========================================

try:
    rmtree(save_directory)
except:
    pass

print("\nStarting simulations using {} CPUs and the following event probabilities:".format(threadNum))
print("  * inversion:     {:5.2f} %".format(100.0*probas["Inversion"]))
print("  * translocation: {:5.2f} %".format(100.0*probas["Translocation"]))
print("  * duplication:   {:5.2f} %".format(100.0*probas["Duplication"]))
print("  * deletion:      {:5.2f} %".format(100.0*probas["Deletion"]))
print("  * fusion:        {:5.2f} %".format(100.0*probas["Fusion"]))
print("  * fission:       {:5.2f} %".format(100.0*probas["Fission"]))
print("  * WGD:           {:5.2f} %".format(100.0*probas["WGD"]))

durations = {}
for simNumber in range(1,nbSimulations+1):
    started = time.time()

    # normalize probabilities
    total = sum(probas.values())
    if total>1:
        for param in probas.keys():
            probas[param] = probas[param]/total

    # print("---------------------- define functions ---------------------")
    # print("[+] define evolution functions")
    # print("\t[*] tree generator")
    from biological_functions import make_tree
    # print("\t[*] inversion")
    from tme3 import inversion
    # print("\t[*] translocation")
    from biological_functions import reciprocal_translocation
    # print("\t[*] duplication")
    from biological_functions import duplication
    # print("\t[*] deletion")
    from tme3 import deletion
    # print("\t[*] fusion")
    from biological_functions import fusion
    # print("\t[*] fission")
    from tme3 import fission
    # print("\t[*] whole genome duplication")
    from biological_functions import wgd
    # print("\t[*] other functions")
    from biological_functions import print_genome, make_root_genome

    #print("\n--------------------- initializing pool ----------------------")
    waiting_queue    = Queue()
    the_results      = Queue()
    poison_pill      = Queue()
    transition_space = Queue()
    l=Lock()

    #print("[+] defining workers on "+str(cpu_count())+" available processor(s)")
    def worker_main(queue):
        while True:
            item = queue.get(True)
            eventCount=item["counts"]
            # get genome data and apply events
            genome=item["data"]
            print(genome)
            for event in item["rearrangements"]:
                if event=="Inversion":
                    genome=inversion(genome, parameters["mean_inv_len"])
                    eventCount["inversions"] += 1
                elif event=="Translocation":
                    genome=reciprocal_translocation(genome)
                    eventCount["translocations"] += 1
                elif event=="Duplication":
                    genome=duplication(genome, parameters["mean_dup_len"])
                    eventCount["duplications"] += 1
                elif event=="Deletion":
                    genome=deletion(genome, parameters["mean_del_len"])
                    eventCount["deletions"] += 1
                elif event=="Fusion":
                    genome=fusion(genome)
                    eventCount["fusions"] += 1
                elif event=="Fission":
                    genome=fission(genome)
                    eventCount["fissions"] += 1
                elif event=="WGD":
                    genome=wgd(genome, parameters["averageDeletionRateWGD"])
                    eventCount["wgds"] += 1
                else:
                    pass
                #print(print_genome(genome))
            item["data"] = genome
            item["counts"] = eventCount
            transition_space.put(item)

    pool = Pool(threadNum, worker_main,(waiting_queue,)) # create pool of worker functions
    # be careful not to forget the comma here -------^

    print("\n--------------------- Simulation {} ----------------------".format(simNumber))

    t = make_tree(nbSpecies)
    t.get_tree_root().name="root"
    #print(t.write(format=1))
    #print(t)

    root=t.get_tree_root()
    root = {}
    root["name"] = t.get_tree_root().name
    root["data"] = make_root_genome(nbChroms, nbGenes)
    root["rearrangements"] = [None]
    root["counts"]={ "inversions":0,
        "translocations":0,
        "duplications":0,
        "deletions":0,
        "fusions":0,
        "fissions":0,
        "wgds":0
    }

    print("root genome has "+str(nbGenes)+" gene(s) on "+str(nbChroms)+" chromosome(s):")
    print("{"),
    tmplist = []
    for ch in sorted(root["data"].keys()):
        tmplist.append( (str(ch),str(len(root["data"][ch])-1)) )
    print(", ".join([ "{}: {}".format(x[0],x[1]) for x in tmplist ])),
    print("}")
    # initialize computation by placing the first task in the queue
    transition_space.put(root)

    # print("\n---------------------- computing evol ----------------------")
    nbTasks = len([node for node in t.traverse()])
    computed=[]
    while(len(computed)<nbTasks):
        # get finished tasks from transition space and launch next tasks

        descendants=[]
        report = transition_space.get(True)
        # msg="[-] Task "+str(report["name"])+" has finished. "

        try:
            descendants = [child.name for child in t.search_nodes(name=report["name"])[0].get_children()]
            # if descendants!=[]:
            #     msg=msg+"Begin children "+", ".join([str(d) for d in descendants])
            # else:
            #     msg=msg+"No descendants (leaf)"

            for descendant in descendants:
                task={}
                task["name"]=descendant
                task["data"]=report["data"]
                task["parameters"]=parameters
                task["counts"]={ "inversions":0,
                    "translocations":0,
                    "duplications":0,
                    "deletions":0,
                    "fusions":0,
                    "fissions":0,
                    "wgds":0
                }

                # generate rearrangements here
                phylo_distance = t.get_distance(report["name"],descendant)
                rea_on_branch= int(E * phylo_distance)
                if rea_on_branch<minimalEventsOnBranch:
                    rea_on_branch=minimalEventsOnBranch
                # then we choose the type of events
                rearrangements=[]
                rearrangements = list(np.random.choice(sorted(probas.keys()), rea_on_branch, p=[probas[event] for event in sorted(probas.keys())]))
                shuffle(rearrangements)
                task["rearrangements"]=rearrangements

                waiting_queue.put_nowait(task)
        except:
            # msg=msg+"No descendants (leaf)"
            pass

        #l.acquire()
        #print(msg)
        #l.release()

        computed.append(report)
        the_results.put(report)

    pool.close() # no more tasks to add


    # wait for last tasks to finish
    results=[]
    while len(results)<nbTasks:
        try:
            results.append(the_results.get())
        except:
            pass


    if not exists(save_directory):
        makedirs(save_directory)
    if not exists(save_directory+"/"+str(simNumber)):
        makedirs(save_directory+"/"+str(simNumber))

    numInversions = 0
    numTranslocations = 0
    numDuplications = 0
    numDeletions = 0
    numFusions = 0
    numFissions = 0
    numWgds = 0

    with open(save_directory+"/"+str(simNumber)+"/genomes.txt","w") as out :
        for r in results:
            numInversions += r["counts"]["inversions"]
            numTranslocations += r["counts"]["translocations"]
            numDuplications += r["counts"]["duplications"]
            numDeletions += r["counts"]["deletions"]
            numFusions += r["counts"]["fusions"]
            numFissions += r["counts"]["fissions"]
            numWgds += r["counts"]["wgds"]
            if len(r["name"]) > 0 and r["name"][0] != 'a':
                out.write(">"+r["name"]+"\n")
                for ch in sorted(r["data"].keys()):
                    out.write(" ".join([str(elt) for elt in r["data"][ch] if elt != 0 ])+" $\n")

    # print("\t-tree file")
    with open(save_directory+"/"+str(simNumber)+"/tree-newick.txt","w") as outtree:
        outtree.write(t.write(format=1))

    duration = str(datetime.timedelta(seconds=time.time()-started))
    durations[simNumber] = duration.replace(":","h ",1).replace(":","m ",1)+"s"

    with open(save_directory+"/"+str(simNumber)+"/statistics.txt","w") as outStat:
        outStat.write("--------------------- Simulation {} ----------------------\n".format(simNumber))
        outStat.write("root genome has "+str(nbGenes)+" gene(s) on "+str(nbChroms)+" chromosome(s):\n")
        outStat.write("{ ")
        tmplist = []
        for ch in sorted(root["data"].keys()):
            tmplist.append( (str(ch),str(len(root["data"][ch])-1)) )
        outStat.write(", ".join([ "{}: {}".format(x[0],x[1]) for x in tmplist ]))
        outStat.write("}\n")
        outStat.write("---------------------- Statistics -----------------------\n")
        outStat.write("Inversions:     {}\n".format(numInversions))
        outStat.write("Translocations: {}\n".format(numTranslocations))
        outStat.write("Duplications:   {}\n".format(numDuplications))
        outStat.write("Deletions:      {}\n".format(numDeletions))
        outStat.write("Fusions:        {}\n".format(numFusions))
        outStat.write("Fissions:       {}\n".format(numFissions))
        outStat.write("WGD:            {}\n".format(numWgds))

    print("---------------------- Statistics -----------------------".format(simNumber))
    print("Inversions:     {}".format(numInversions))
    print("Translocations: {}".format(numTranslocations))
    print("Duplications:   {}".format(numDuplications))
    print("Deletions:      {}".format(numDeletions))
    print("Fusions:        {}".format(numFusions))
    print("Fissions:       {}".format(numFissions))
    print("WGD:            {}".format(numWgds))
    print("") #print("---------------------------------------------------------")
    print("Time: {}\n".format(durations[simNumber]))

#======================================== main simulation loop ========================================
# for simNumber in sorted(durations.keys()):
#     print("simulation"),
#     print(simNumber),
#     print(":"),
#     print(durations[simNumber])
