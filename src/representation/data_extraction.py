import csv 
import matplotlib.pyplot as plt
import dadi
import math
import numpy as np
import argparse
import time 


# ---------------------------------------
# ---------- Goal of this file ----------
# ---------------------------------------
# The goal of this file is to extract the dadi and stairway results we have
# from the various output file format and save them in csv files 

# ---------------------------------------


# -----------------------------------------------
# ------------------ Arguments ------------------
# -----------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--dadi_louis', type=str, required=True)
parser.add_argument('--stairway_path', type=str, required=True)

parser.add_argument('--nbr_generation', type=int, required=True) # number of generation used in the simulation
parser.add_argument('--Ninit', type=int, required=True) # Initial popsize

parser.add_argument('--simu_type', type=str, required=True)  # simulation type (string like expansion_neutral)

parser.add_argument('--rep', type=int, required=True)  # simulation type (string like expansion_neutral)
parser.add_argument('--GR', type=int, required=True)  # simulation type (string like expansion_neutral)


args = parser.parse_args()

dadi_louis_path = args.dadi_louis # Path to the csv file containning the results of Louis' dadi script
stairway_dfh_path = args.stairway_path # Path to the Stairway output file of DFH
type_simu = args.simu_type
GR = args.GR
rep = args.rep

directory_name = ['../data/dadi/coords/'+ str(rep),'../data/stairway/coords/'+ str(rep), '../data/dadi/'+ str(rep), '../data/stairway/'+ str(rep) ]

for dname in directory_name :
	try:
		os.mkdir(dname)
		print(f"Directory '{dname}' created successfully.")
	except FileExistsError:
		print(f"Directory '{dname}' already exists.")



# -------------------------------------------

# Theorical values
nbr_gen = args.nbr_generation
Popsize_init = args.Ninit


# ----------------------------------------------------------------------------------
# -------------------------------------- DADI --------------------------------------
# ----------------------------------------------------------------------------------

Nu_dadi = []
T_dadi = []
# We gather the values of dadi from Louis script
csv_file = open(dadi_louis_path)
csv_content = csv.reader(csv_file, delimiter = ",")

for ligne in csv_content :
	Nu_dadi.append(float(ligne[0]))
	T_dadi.append(float(ligne[1]))

# We calculate the MEDIAN (and not the mean to lower the impact of anomalies) of the infered params
mean_Nu_dadi = np.median(Nu_dadi)
mean_T_dadi = np.median(T_dadi)

# List of lists containing the coordinates of the points
data_dadi = [ ['x', 'y'] ]
# The generations
time_gen = [x for x in range(0, nbr_gen + 5, 5)]

for t in time_gen :  # for every generation
	
	if t >= (mean_T_dadi  * (2*Popsize_init)) : # if it is before the event the population size has not changed 
		data_dadi.append([float(t) / (2*Popsize_init), Popsize_init]) 
	else : # if it is after the event, the population size has changed depending on the simulation type
		if 'reduction' in type_simu :
			data_dadi.append([float(t) / (2*Popsize_init), Popsize_init * mean_Nu_dadi])
		else : 
			if 'expansion' in type_simu :# 'expansion'
				data_dadi.append([float(t) / (2*Popsize_init) , Popsize_init * np.exp(np.log(mean_Nu_dadi)/mean_T_dadi*(mean_T_dadi - (t / (2*Popsize_init)) ) )])
			else : # constant
				data_dadi([float(t) / (2*Popsize_init), Popsize_init])
# -------------------------------------------------------
# ------------------ CSV FILE CREATION ------------------
# -------------------------------------------------------

# The data, the coordinates (x : time in generation, y : popsize)
with open('../data/dadi/coords/'+ str(rep) + "/" + type_simu + "_" + str(GR) + '.csv', 'w', newline='') as csvfile:
	writer = csv.writer(csvfile)
	writer.writerows(data_dadi)
# Csv inferred params (T, Nu)
inf_params = [ ['T', 'Nu'], [mean_T_dadi, mean_Nu_dadi] ]
with open('../data/dadi/' + str(rep) + "/infered_params_dadi_" + type_simu + "_" + str(GR) + '.csv', 'w', newline='') as csvfile:
	writer = csv.writer(csvfile)
	writer.writerows(inf_params)


# ----------------------------------------------------------------------------------
# ------------------------------------ STAIRWAY ------------------------------------
# ----------------------------------------------------------------------------------

# Step 1 : We gather and sort the data
# The Nu estimated for every sexual reproduction frequency
Nu_stairway = []
# Popsize for a few generations after the event (at the end)
Ne_fin = []
# Popsize for a few generations before the event (at the start)
Ne_init = []
# T estimated (with the growth or decrease of the popsize slope)
T_stairway = []

# List of lists containing the coordinates of the points
data = []
dic_t_ne = {} # t keys, ne values
with open(stairway_dfh_path) as input_file:
	line = input_file.readline() # we skip the first line (the columns' names)
	line = input_file.readline()
	while line != '':
		
		Ne = float(line.split('\t')[6])
		t = float(line.split('\t')[5])
		if t not in list(dic_t_ne.keys()) :  # If we don't already have a y for that x
			data.append([t, Ne])
			#print(t,Ne)
			dic_t_ne[t] = Ne 
		line = input_file.readline()
data.sort()
data.insert(0, ['x', 'y'])

# -------------------------------------------------------
# ------------------ CSV FILE CREATION -- 1st step ------
# -------------------------------------------------------
# The data, the coordinates (x : time, y : popsize)
with open('../data/stairway/coords/' + str(rep) + '/data_stairway_' + type_simu + "_" + str(GR) + '.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data)
# -------------------------------------------------------

# -------------------------------------------------------
# ------------------ DETERMINE T & Nu -------------------
# -------------------------------------------------------
list_t_ne = sorted(dic_t_ne.items(), reverse=True) # we sort the elements (the time) from the oldest to the newest

# We calculate T
for error_margin in range(1, 10) : # for various error margin
	i = 0 # starting point (the oldest time)
	interval_erreur_accept = [error_margin/1000, 1 + error_margin/1000] # accepted variation (growth or decrease)
	while i + 1 < len(list_t_ne) : # we go through all the points (from oldest to newest)
		(tv, ne_v) = list_t_ne[i]
		(tj, ne_j) = list_t_ne[i+1]
		rapport = float(ne_v) / float(ne_j)
		
		if rapport >= interval_erreur_accept[0] and rapport <= interval_erreur_accept[1] : # in the accepeted error evolution interval
			i = i + 1
		else :
			rapports = [rapport]
			for j in range(2, 7) : # we test if it is not an anomaly ( if the trend is durable )
				(tjj, ne_jj) = list_t_ne[i+j]
				rapports.append(float(ne_v) / float(ne_jj))
			mean_rap = np.mean(rapports)
			if mean_rap >= interval_erreur_accept[0] and mean_rap <= interval_erreur_accept[1] : # in the accepeted error evolution interval
				i = i + 1
			else : # if it wasn't an anomality, we keep the t
				T_stairway.append(tv)
				i = len(list_t_ne) # we get out of the while loop


# We define the Infered stairway T as the mean of the various possibles T we had
T_mean = np.mean(T_stairway)

# We calculate Nu

for t, ne in list_t_ne[0:10] : # We go through the 10 first generations (to get the average Ninit)
	if t <= T_mean : # after the event
		Ne_fin.append(ne)
	else :
		Ne_init.append(ne)
      
for t, ne in list_t_ne[int(len(list_t_ne)-20): int(len(list_t_ne)-10)] : # We go through the 10 to 20 last generations (to get the average Nfinal)
	if t <= T_mean : # after the event
		Ne_fin.append(ne)
	else :
		Ne_init.append(ne)


def moyenne_geometrique(liste):
    """ Returns the geometric mean of a given list """
    if not liste:
        return None
    if any(x <= 0 for x in liste):
        raise ValueError("Tous les éléments doivent être strictement positifs.")
    produit = math.prod(liste)  
    return produit ** (1 / len(liste))

# We determine the mean Ninit and Nfinal
Ne_init_mean = moyenne_geometrique(Ne_init)
Ne_fin_mean = moyenne_geometrique(Ne_fin)

print('Mean ----------------------------------------------------- \n', Ne_init, Ne_fin)

# We calculate the infered Nu stairway as the ratio of Nfinal / Ninit
Nu_stairway = Ne_fin_mean / Ne_init_mean


# -------------------------------------------------------
# ------------------ CSV FILE CREATION -- 2nd step ------
# -------------------------------------------------------
# The infered parameters T, Nu
inf_params = [ ['T', 'Nu'], [T_mean, Nu_stairway] ]
with open('../data/stairway/' + str(rep) + '/infered_params_stairway_' + type_simu + "_" + str(GR) + '.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(inf_params)