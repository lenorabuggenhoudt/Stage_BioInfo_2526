#!/usr/bin/python

## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 

# This program will use the .trees files (= tree sequence, TS) given by SLiM to compute summary statistics (pi)

import msprime, pyslim, tskit, numpy as np, sys, warnings
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--replicatenumber', type=int, required=True)
parser.add_argument('--nb_samples', type=int, required=True)
parser.add_argument('--mutation_rate', type=float, required=True)
parser.add_argument('--popsize', type=int, required=True)
parser.add_argument('--windows', type=int, required=True)
parser.add_argument('--recombinationrate', type=float, required=True)
parser.add_argument('--generationmitoses', type=int, required=True)
parser.add_argument('--out_dir', type=str, required=True)
args = parser.parse_args()


replicatenumber = args.replicatenumber # Replicate number
nb_samples = args.nb_samples # Number of sample for the simplication and sumstats (diploïd) 
popsize = args.popsize # Size of the population 
mutation_rate = args.mutation_rate # Mutation rate
recombinationrate = args.recombinationrate # Recombination rate
windows = args.windows
generationmitoses = args.generationmitoses 
out_dir = args.out_dir



##### SAMPLING, ADDING NEUTRAL MUTATION AND GETTING SNP MATRIX #####
##### RECAPITATE, SAMPLING AND ADDING NEUTRAL MUTATION) ###################################
# https://tskit.dev/pyslim/docs/latest/tutorial.html

ts = tskit.load(f"./{out_dir}/ts_{replicatenumber}_{generationmitoses}.trees")
ts = pyslim.update(ts)

warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

# Recapitate the simulation to provide a “prior history” for the initial generation of the simulation = fully coalesced #
# RTS : recapitated 

rts = pyslim.recapitate(ts, recombination_rate=recombinationrate*generationmitoses, ancestral_Ne=popsize)
      
# Sampling indivuals from the TS to compute sumstats later # 
# SRTS = simplified recapitated TS 

rng = np.random.default_rng()
alive_inds = pyslim.individuals_alive_at(rts, 0) # sample individuals alive = coming from the last generation
if len(alive_inds) < nb_samples : # If we have less individuals alive than number of samples asked for
  warnings.warn('The number of samples the programm was asked to keep is greater than the final size of the population (after simulation), the programm will use all the samples of the final population \n')
  nb_samples = len(alive_inds) # we will use all of the individuals at our disposition 
keep_indivs = rng.choice(alive_inds, nb_samples, replace=False)
keep_nodes = []
for i in keep_indivs:
  keep_nodes.extend(rts.individual(i).nodes)

# sort the list ???                                                                   ############################# /!\ /!\ /!\ /!\
list.sort(keep_nodes)


srts = rts.simplify(keep_nodes, keep_input_roots=True)

# Adding neutral mutations to the samples (overlay) # 
# OSRTS = overlay simplified recapitated TS 

osrts = msprime.mutate(srts, rate=(mutation_rate), keep=True) # keep existing mutations

##### SUMMARY STATISTICS (WITH WINDOWS) ####################################

# We want a user defined number of windows for the first chromosome but only 1 value
# for the second since it's entirely neutral 

windows_chr = np.linspace(0, osrts.sequence_length/2, num=windows)
windows_chr = [int(x) for x in windows_chr] + [int(osrts.sequence_length)]

## PI ##
pi = osrts.diversity(windows=windows_chr, mode='site')
pi_chr1 = pi[:-1]
pi_chr2 = pi[-1]
#pi_chr1.tofile(f"./{out_dir}/pi_chr1_{replicatenumber}_{generationmitoses}.csv", sep=',')
#pi_chr2.tofile(f"./{out_dir}/pi_chr2_{replicatenumber}_{generationmitoses}.csv", sep=',')


# add other stats here ?
# does not work treesequence not iterable
"""
robinson_foulds = [0] # test to see if something is written
for tree in osrts :
    ligne = []
    for tree2 in osrts :
      rf = tree.rf_distance(tree2)
      ligne.append(rf)
    robinson_foulds.append(ligne)
robinson_foulds.tofile(f"results/test_robinson_foulds.csv", sep=',')
"""

# ld
ld = tskit.LdCalculator(osrts)
r2 = ld.r2_matrix()
#r2.tofile(f"./{out_dir}/test_ld_r2_matrix_{replicatenumber}_{generationmitoses}.csv", sep =',')


#sfs
site_freq = osrts.allele_frequency_spectrum(windows=windows_chr)
# windows ??
#site_freq_1 = site_freq[:-1]
#site_freq_2 = site_freq[-1]
#site_freq_1.tofile(f"./{out_dir}/test_sfs_1_{replicatenumber}_{generationmitoses}.csv", sep=',')
#site_freq_2.tofile(f"./{out_dir}/test_sfs_2_{replicatenumber}_{generationmitoses}.csv", sep=',')



print(osrts.draw_text())