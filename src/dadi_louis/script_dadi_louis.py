#!/usr/bin/python

## Author : Louis OLLIVIER ~ louis.ollivier@etu.univ-rouen.fr , modified by Lénora BUGGENHOUDT	~ lenora.buggenhoudt@universite-paris-saclay.fr
## Université Paris-Saclay
## Lab : LISN ~ UMR9015 ~ BIOINFO team 


# and then do a demography parameters inference with dadi (thanks to the SFS).

import numpy as np, dadi, random
import argparse
import time # for the timestamp 
import os

parser = argparse.ArgumentParser()
parser.add_argument('--config_file', type=str, required=True)
parser.add_argument('--GR', type=int, required=True)
parser.add_argument('--repetiton_number', type=int, required=True)
#parser.add_argument('--output', type=str, required=True)
args = parser.parse_args()

# -------------------------------------------------
# -------------------------------------------------
import yaml
with open(args.config_file, 'r') as f:
    data = yaml.full_load(f)
sfs_file = data.get('path_to_sfs')

if 'expansion' in data.get('vcf') :
	simu_type = 'expansion_'
	lower_bound = [1, 0.001]
	upper_bound = [50, 1]
else :
	simu_type = 'reduction_'
	lower_bound = [0.001, 0.001] 
	upper_bound = [1, 1]
	
if 'neutral' in data.get('vcf') :
    simu_type = simu_type + 'neutral'
else :
    simu_type = simu_type + 'sweep'
# --------------------------------------------------

directory_name = simu_type + "_" + str(args.repetiton_number)

# Create the directory
try:
    os.mkdir("../../results/dadi_louis/" + directory_name )
    print(f"Directory '{directory_name}' created successfully.")
except FileExistsError:
    print(f"Directory '{directory_name}' already exists.")


# ------------------------------------------------------------------
##### DEMOGRAPHY INFERENCE #####
ns=5  #n
pts = [40] # § Implementation : https://dadi.readthedocs.io/en/latest/user-guide/simulation-and-fitting/
# The number of grid points in each dimension for representing ϕ

f = "../deminfhelper/" + sfs_file
data = dadi.Spectrum.from_file(f)#dadi.Spectrum(rawSFS) 

func = dadi.Demographics1D.growth
func_ex = dadi.Numerics.make_extrap_log_func(func) # Make the extrapolating version of our demographic model function.

# IIntervals for the parameters to be infered

# True values : Nu = Ncontemp/Ninit = 15,747 ; T = 1,5 (2Na)
#lower_bound = [0.01, 0.01] 
#upper_bound = [50, 50]  

inf_param = [] # Will contain the values for the infered parameters 

# ~ thetaW_dadi = data.Watterson_theta()
# ~ print("Wtheta dadi :", thetaW_dadi)

for j in range(4) :
	for i in range(0,5): # We will do the parameters inference many times
			
		# Initial guess for the parameters, because we know the true parameters will pick them randomly in the intervals
		p0 = [random.uniform(lower_bound[0],upper_bound[0]),
			random.uniform(lower_bound[1],upper_bound[1])]

		# Perturb our parameters before optimization. This does so by taking each parameter a up to a factor of two up or down.
		p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound) 

		# ~ print("Beginning optimization ************************************************")
		popt = dadi.Inference.optimize_log_fmin(p0, data, func_ex, pts, 
														lower_bound=lower_bound, 
														upper_bound=upper_bound,
														maxiter=50) # verbose=len(p0),
		# ~ print("Finished optimization **************************************************")
		
		inf_param.append(popt)
	
print(inf_param)
#np.savetxt("infered_param_"  + str(time.time()) + ".csv", inf_param, delimiter=',')
np.savetxt("../../results/dadi_louis/" + directory_name + "/infered_param_" + simu_type + "_" + str(args.GR) + "_" + ".csv", inf_param, delimiter=',')

# ~ https://dadi.readthedocs.io/en/latest/examples/YRI_CEU/YRI_CEU/
# ~ https://dadi.readthedocs.io/en/latest/user-guide/specifying-a-model/	
# ~ §§ voir demo graphy inference suite 
