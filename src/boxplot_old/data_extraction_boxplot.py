import csv 
import matplotlib.pyplot as plt
import dadi
import math
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

print(Nu_dadi)
print(T_dadi)

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
Nu_stairway = []
Nu_fin = []
Nu_init = []
T_stairway = []
with open(stairway_dfh_file) as input_file:
	line = input_file.readline() # we skip the first line (the columns' names)
	line = input_file.readline()
	while line != '':
		nu = float(line.split('\t')[6])
		t = float(line.split('\t')[5])
		if t <= T_theo :
			Nu_fin.append(nu)
		else :
			Nu_init.append(nu)
		Nu_stairway.append(nu)
		T_stairway.append(t)
		line = input_file.readline()
      
popfinale = 10000
Nu_stairway_rap = [a / b for a, b in zip(Nu_fin, Nu_init)]
print('Rap -----------------------------------------------------')
#print(Nu_stairway_rap)
Nu_stairway_geo_mean = []
"""
for i in range(0, len(Nu_stairway_rap)) :
	for j in range(i+1,len(Nu_stairway_rap) ) :
		gm = math.sqrt(Nu_stairway_rap[i] * Nu_stairway_rap[j])
		Nu_stairway_geo_mean.append(gm)
"""
i = 0
while i + 1 < (len(Nu_stairway_rap)) :
	gm = math.sqrt(Nu_stairway_rap[i] * Nu_stairway_rap[i + 1])
	Nu_stairway_geo_mean.append(gm)
	i = i + 2
print('Mean ----------------------------------------------------- \n', Nu_stairway_geo_mean)
Nu_stairway = Nu_stairway_geo_mean
#Nu_stairway = [ popfinale/x for x in Nu_stairway]  
            
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
plt.savefig("../../results/boxplots/" + "boxplot_nu_" + type_simu + "_" + str(time.time()) + ".pdf")              
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
plt.savefig("../../results/boxplots/" + "boxplot_t_" + type_simu + "_" + str(time.time()) + ".pdf")                 
#plt.show()


