# Themis' script modified
import tskit 
import sys
import os
import numpy as np
import msprime, pyslim, warnings
import argparse
import time

parser = argparse.ArgumentParser()
parser.add_argument('--replicatenumber', type=int, required=True)
parser.add_argument('--nb_samples', type=int, required=True)
parser.add_argument('--mutation_rate', type=float, required=True)
parser.add_argument('--popsize', type=int, required=True)
parser.add_argument('--recombinationrate', type=float, required=True)
parser.add_argument('--generationmitoses', type=int, required=True)
parser.add_argument('--out_dir', type=str, required=True)
parser.add_argument('--reductionratio', type=int, required=True)
parser.add_argument('--expansionratio', type=int, required=True)
args = parser.parse_args()


replicatenumber = args.replicatenumber # Replicate number
nb_samples = args.nb_samples # Number of sample for the simplication and sumstats (diploïd) 
popsize = args.popsize # Size of the population 
mutation_rate = args.mutation_rate # Mutation rate
recombinationrate = args.recombinationrate # Recombination rate
generationmitoses = args.generationmitoses 
reductionratio = args.reductionratio
expansionratio = args.expansionratio
out_dir = args.out_dir

# Load the trees
full_tree = tskit.load(f"./{out_dir}/ts_{replicatenumber}_{generationmitoses}.trees") # to change, also need to not supress the tree file before calling this file
warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)
nbr_sites_pre = full_tree.num_sites
print(f"Total number of variant sites pre adding neutral mutations : {nbr_sites_pre}")


"""
# Create individual names based on their population
individual_names = []
for individual in full_tree.individuals():
    node_id = individual.nodes[0]  # Assuming 1 node per individual
    pop_id = full_tree.node(node_id).population
    if pop_id == 0:
        individual_names.append(f"p1_ind{individual.id}")
    elif pop_id == 1:
        individual_names.append(f"p2_ind{individual.id}")
    else:
        individual_names.append(f"unknown_individual{individual.id}")
"""
# Simplify the tree sequence with the samples of interest
#simplified_ts = full_tree.simplify(samples=full_tree.samples())

# Recapitate the simulation to provide a “prior history” for the initial generation of the simulation = fully coalesced #
# RTS : recapitated 
rts = pyslim.recapitate(full_tree, recombination_rate=recombinationrate*generationmitoses, ancestral_Ne=popsize)

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



# SRTS = simplified recapitated TS 
srts = rts.simplify(keep_nodes, keep_input_roots=True)
# Adding neutral mutations to the samples (overlay) # 
# OSRTS = overlay simplified recapitated TS 

#osrts = msprime.mutate(srts, rate=(mutation_rate), keep=True) # keep existing mutations
osrts = msprime.sim_mutations(srts, rate=(mutation_rate), keep=True) # keep existing mutations

nbr_sites_post = osrts.num_sites
print(f"Total number of variant sites post adding neutral mutations : {nbr_sites_post}")


# Site frequency spectrum
#results = osrts.allele_frequency_spectrum(sample_sets=None, windows=None, mode='site', span_normalise=False, polarised=False)
results = osrts.allele_frequency_spectrum(sample_sets=None, windows=None, mode='site', polarised=False)
normalized_results = results / results.sum()
print(normalized_results.sum())  # ça donnera bien 1.0
results = normalized_results
# polarised False --> folded SFS
results[0] = nbr_sites_post # the first element is not interesting for us (it does not correspond to the singleton, so we replace it by the number of sites
results.tofile(f"./{out_dir}/sfs_{replicatenumber}_{generationmitoses}.csv", sep=',')  
# Get vcf from tree
#os.makedirs(f"./{out_dir}/"+"/vcf/"+f"GR{generationmitoses}/",exist_ok=True) # if we want the GR files
time_indicator = time.time()
os.makedirs(f"./{out_dir}/"+"/vcf/",exist_ok=True)
#output_gr_directories =f"./{out_dir}/"+"/vcf/" +f"GR{generationmitoses}/" + f"ts_{replicatenumber}_{generationmitoses}_{time_indicator}" + ".vcf"
output_timed =f"./{out_dir}/"+"/vcf/"+ f"ts_{replicatenumber}_{generationmitoses}_{time_indicator}" + ".vcf"
output=f"./{out_dir}/"+"/vcf/"+ f"ts_{replicatenumber}_{generationmitoses}" + ".vcf"
#summary_file_gr_directories = f"./{out_dir}/"+"/vcf/" +f"GR{generationmitoses}/" + "summary_file.txt"
summary_file = f"./{out_dir}/"+"/vcf/"  + f"summary_file_{replicatenumber}_{generationmitoses}.txt"
summary_file_timed = f"./{out_dir}/"+"/vcf/"  + f"summary_file_{replicatenumber}_{generationmitoses}_{time_indicator}.txt"
with open(output,"w") as file:
    #simplified_ts.write_vcf(file, individual_names=individual_names,position_transform=lambda x: np.fmax(1,x))
    osrts.write_vcf(file,position_transform=lambda x: np.fmax(1,x)) # bc we only have one population, the individual names based on the population don't matter 

#with open(output_timed,"w") as file:
    #simplified_ts.write_vcf(file, individual_names=individual_names,position_transform=lambda x: np.fmax(1,x))
    #osrts.write_vcf(file,position_transform=lambda x: np.fmax(1,x)) # bc we only have one population, the individual names based on the population don't matter 


with open(summary_file, "w") as sum_file:
    sum_file.write("Pop size initiale : " + str(popsize) + "\n")
    sum_file.write("Pop size finale : " + str(len(alive_inds)) + "\n")
    sum_file.write("Number of samples kept : " + str(nb_samples) + "\n")
    sum_file.write("Mutation rate : " + str(mutation_rate) + "\n")
    sum_file.write("Number of sites mutated pre added neutral mutations : " + str(nbr_sites_pre) + "\n")
    sum_file.write("Number of sites mutated post added neutral mutations : " + str(nbr_sites_post) + "\n")
    sum_file.write("Recombination rate : " + str(recombinationrate) + "\n")
    sum_file.write("Generation mitoses (1/alpha, GR) : " + str(generationmitoses) + "\n")

    if "reduction" in out_dir:
        sum_file.write("Type : Reduction (lasts 2000 epochs, start at 8000 epochs), ratio : 1/" + str(reductionratio) +"\n")
        
    else:
            sum_file.write("Type : Expansion (lasts 2000 epochs, start at 8000 epochs), ratio : " + str(expansionratio) +"\n")

    if "neutral" in out_dir:
        sum_file.write("Selection : None\n")
    else:
        sum_file.write("Selection : Yes (start at 2000 epochs)\n")
    sum_file.close()
    