import csv 
import matplotlib.pyplot as plt
import dadi
import math
import numpy as np
import argparse
import time 
# -------------------------------------------
# ------------------ Arguments --------------
# -------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--dadi_louis', type=str, required=True)
#parser.add_argument('--dadi_dfh', type=str, required=True)
parser.add_argument('--stairway_summary', type=str, required=True)
parser.add_argument('--Nu_theo', type=float, required=True)
parser.add_argument('--T_theo', type=float, required=True)
parser.add_argument('--simu_type', type=str, required=True)
args = parser.parse_args()

dadi_louis_file = args.dadi_louis # Path to the csv file containning the results of Louis' dadi script
#dadi_dfh_file = args.dadi_dfh # Path to the dadi output file of DFH
stairway_dfh_file = args.stairway_summary # Path to the Stairway output file of DFH
type_simu = args.simu_type
# -------------------------------------------

# Theorical values (put randomly for now)
Nu_theo = args.Nu_theo
T_theo = args.T_theo

# -------------- dadi values --------------------
# dadi Louis
Nu_dadi = []
T_dadi = []
# dadi DFH
Nu_dadi_dfh = []
T_dadi_dfh = [] 

# We gather the values of dadi from Louis script
csv_file = open(dadi_louis_file)
csv_content = csv.reader(csv_file, delimiter = ",")
print(csv_content)
for ligne in csv_content :
    print(ligne)
    Nu_dadi.append(float(ligne[0]))
    T_dadi.append(float(ligne[1]))
print("Dadi----------")
print("Nu : ", Nu_dadi)
print("T : ", T_dadi)
mean_Nu_dadi = np.mean(Nu_dadi)
mean_T_dadi = np.mean(T_dadi)
data_dadi = [ ['x', 'y'] ]
nbr_gen = 10000
Popsize_init = 10000
time_gen = np.linspace(0, nbr_gen)
for t in time_gen :
	if t >= mean_T_dadi :
		data_dadi.append([t / (2*Popsize_init), Popsize_init])
	else : 
		data_dadi.append([t / (2*Popsize_init), Popsize_init * mean_Nu_dadi])

# ------------------ csv file creation ------------------
with open('../data/dadi/coords/data_dadi.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data_dadi)
# --------------------- csv inferred params ------------------------------
inf_params = [ ['T', 'Nu'], [mean_T_dadi, mean_Nu_dadi] ]
with open('../data/dadi/infered_params_dadi.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(inf_params)
# ---------------------------------------------------------
# We gather the values of dadi from the DFH output file

#dadi_vals = dadi_output_parse('1.InferDM.bestfits')
"""
dadi_file = open(dadi_dfh_file)
#for elem in dadi_vals:
		# elem is structured like this:
		# [iter_number, log_likelihood, [nuB, nuF, TB, TF], theta]
	# Nu_dadi2.append(elem[2][0])
# TRES TRES BROUILLON
for ligne in dadi_file :
	if ligne[0] !='#' and ligne[0] != ' ':
		val = ligne.split()
		if len(val) > 0 :
			Nu_dadi_dfh.append(float(val[1]))
			T_dadi_dfh.append(float(val[2]))
"""			


# -------------------- Stairway values --------------------------- (dfh)
# Step 1 : We gather and sort the data
Nu_stairway = []
Ne_fin = []
Ne_init = []
T_stairway = []
dic_t_ne = {}  # key : t, value : Ne_med
data = [ ['x', 'y'] ]
with open(stairway_dfh_file) as input_file:
	line = input_file.readline() # we skip the first line (the columns' names)
	line = input_file.readline()
	while line != '':
		Ne = float(line.split('\t')[6])
		t = float(line.split('\t')[5])
		if t not in list(dic_t_ne.keys()) :
			data.append([t, Ne])
		dic_t_ne[t] = Ne 
		line = input_file.readline()

# ------------------ csv file creation ------------------
with open('../data/stairway/coords/data_stairway.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(data)
# ---------------------------------------------------------


list_t_ne = sorted(dic_t_ne.items(), reverse=True) # we sort the elements from the oldest to the newest
# T --------------
# Step 2 : We compare the data

for error_margin in range(1, 10) :
	i = 0
	interval_erreur_accept = [error_margin/1000, 1 + error_margin/1000] # changement de 5 accepetés (en baisse ou en hausse)
	while i + 1 < len(list_t_ne) :
		(tv, ne_v) = list_t_ne[i]
		(tj, ne_j) = list_t_ne[i+1]
		rapport = float(ne_v) / float(ne_j)
		
		if rapport >= interval_erreur_accept[0] and rapport <= interval_erreur_accept[1] : # in the accepeted error ecolution interval
			i = i + 1
		else :
			print("Tv : ", tv)
			print("Tj : ", tj)
			print("Rapport : ", rapport)
			rapports = [rapport]
			for j in range(2, 7) : # we test if it is not an anomalie ( if the trend is durable )
				(tjj, ne_jj) = list_t_ne[i+j]
				rapports.append(float(ne_v) / float(ne_jj))
			mean_rap = np.mean(rapports)
			if mean_rap >= interval_erreur_accept[0] and mean_rap <= interval_erreur_accept[1] : # in the accepeted error ecolution interval
				i = i + 1
			else : # if it wasn't an anomality, we keep the t
				T_stairway.append(tv)
				i = len(list_t_ne) # we get out of the while loop

print("\n\nStairway ---------------")
print("T : ", T_stairway)
T_mean = np.mean(T_stairway)
# --------- Nu
"""
for t, nu in list_t_nu :
	if t <= T_mean :
		Nu_fin.append(nu)
	else :
		Nu_init.append(nu)
"""
for t, ne in list_t_ne[0:10] :
	if t <= T_mean :
		Ne_fin.append(ne)
	else :
		Ne_init.append(ne)
      
for t, ne in list_t_ne[int(len(list_t_ne)-20): int(len(list_t_ne)-10)] :
	if t <= T_mean :
		Ne_fin.append(ne)
	else :
		Ne_init.append(ne)


def moyenne_geometrique(liste):
    if not liste:
        return None
    if any(x <= 0 for x in liste):
        raise ValueError("Tous les éléments doivent être strictement positifs.")
    produit = math.prod(liste)  
    return produit ** (1 / len(liste))

Ne_init_mean = moyenne_geometrique(Ne_init)
Ne_fin_mean = moyenne_geometrique(Ne_fin)

print('Mean ----------------------------------------------------- \n', Nu_stairway)
Nu_stairway = Ne_fin_mean / Ne_init_mean
#Nu_stairway = [ popfinale/x for x in Nu_stairway]  

# --------------------- csv inferred params ------------------------------
inf_params = [ ['T', 'Nu'], [T_mean, Nu_stairway] ]
with open('../data/stairway/infered_params_stairway.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(inf_params)

# -------------------------- PLOT ---------------------

# ---------------------- Nu values ----------------------
fig, ax = plt.subplots()
ax.set_ylabel('Nu value')

bplot = ax.boxplot([Nu_dadi, Nu_dadi_dfh, Nu_stairway],
                   patch_artist=True,  # fill with color
                   labels=['Nu_dadi_louis', 'Nu_dadi_dfh', 'Nu_stairway'])  # labels (version python < 3.9), tick_labels (version Python >= 3.9)
      
      
plt.hlines(y=[Nu_theo], xmin=0, xmax=2.5, colors=['r'], linestyles=['--'])             
   


plt.title('Nu param infered')
plt.tight_layout()

#plt.savefig("boxplot_nu_" + type_simu + "_" + str(time.time()) + ".pdf")      
plt.savefig("../../results/boxplots_old/" + "boxplot_nu_" + type_simu + "_" + str(time.time()) + ".pdf")              
#plt.show()



# ------------------------ T values -------------------

fig, ax = plt.subplots()
ax.set_ylabel('T value')


bplot = ax.boxplot([T_dadi, T_dadi_dfh, T_stairway],
                   patch_artist=True,  # fill with color
                   labels=['T_dadi_louis', 'T_dadi_dfh', 'T_stairway'])  # labels (version python < 3.9), tick_labels (version Python >= 3.9)
      
      
plt.hlines(y=[T_theo], xmin=0, xmax=2.5, colors=['r'], linestyles=['--'])             
   


plt.title('T param infered')
plt.tight_layout()
#plt.savefig("boxplot_t_" + type_simu + "_" + str(time.time()) + ".pdf") 
plt.savefig("../../results/boxplots_old/" + "boxplot_t_" + type_simu + "_" + str(time.time()) + ".pdf")                 
#plt.show()



