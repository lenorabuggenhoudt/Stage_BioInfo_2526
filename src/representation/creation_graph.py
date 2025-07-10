import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import argparse
import time 
# -------------------------------------------
# ------------------ Arguments --------------
# -------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--path_dadi_inf_params', type=str, required=True)
parser.add_argument('--path_stair_inf_params', type=str, required=True)
parser.add_argument('--Nu_theo', type=float, required=True)
parser.add_argument('--T_theo', type=float, required=True)
parser.add_argument('--Ninit', type=int, required=True)
parser.add_argument('--simu_type', type=str, required=True)
args = parser.parse_args()

# Inferred params recuperation
df = pd.read_csv(args.path_stair_inf_params) 
T_stair = float(df['T'][0])  # The infered parameters for stairway
Nu_stair = float(df['Nu'][0])

df = pd.read_csv(args.path_dadi_inf_params) 
T_dadi = float(df['T'][0])  # The infered parameters for dadi
Nu_dadi = float(df['Nu'][0])


Nu_theo = args.Nu_theo
T_theo = args.T_theo
Ninit = args.Ninit
type_simu = args.simu_type




# We define a function that given some parameters returns the theorical value of the Ne initial
# generations ---> kind of list of all the generations
# ne ---> taille de popu initiale
# T ---> moment où le changement (expansion ou reduction) a commencé ???? Voir schéma
# nu ----> rapport popu finale / popu initiale



def attendu_expo(generations, Ne, T, nu):
    res = []
    beta = np.log(nu) / T
    
    for gen in generations: # for every generation
        if gen >= T: # if we are before the demographic change, the size of the popu is the same as the initial size (Ne)
            res.append(Ne)
        else:  # if we are during or after the demographic change, we need to calculate the change
            Ne_calc = Ne * np.exp(-beta * (gen - T))
            res.append(Ne_calc)
    
    return np.array(res)


def attendu_reduc(generations, Ne, T, nu) :
    res = []
    beta = np.log(nu) / T
    
    for gen in generations: # for every generation
        if gen >= T: # if we are before the demographic change, the size of the popu is the same as the initial size (Ne)
            res.append(Ne)
        else:  # if we are during or after the demographic change, we need to calculate the change
            Ne_calc = Ne * nu
            res.append(Ne_calc)
    
    return np.array(res)



def attendu(type_simu, generations, Ne, T, nu) :
    if 'reduction' in type_simu :
        return attendu_reduc(generations, Ne, T, nu)
    else :
        return attendu_expo(generations, Ne, T, nu)


# ----------------------------------------------------------------------------

mu = 1e-7
#generations = np.arange(0, 10000, 100) # generations ?
# Adjust generations to cover the full time span (0 to > T_EGR)
T_max = max(T_theo * 2 * 10000, T_dadi * 2 * 10000, T_stair * 2 * 10000)
generations = np.arange(0, int(T_max) + 1000, 100)
T = T_theo * 2 * 10000

# --- Théoriques connues ---
theo_values = attendu(type_simu, generations, Ninit, T, Nu_theo)

theo = pd.DataFrame({
    'model': 'theorique',
    'generations': generations,
    'value': theo_values
})

# --- Dadi ---
# Valeurs estimées et déduite avec dadi # 
# ATTENTION: valeurs assez différentes dans le cas de 1000GR si on prend moyenne ou médiane 
# On prendra la médiane ici
Nfinale = Nu_dadi * Ninit
theta = 4 * Nfinale * mu # 4Nfinale * mu = 4 * (Ninit * ratio_expo) * mu
Ne_dadi = theta / (4 * mu * Nu_dadi) # N initiale
T = T_dadi * 2 * 10000
dadi_values = attendu(type_simu, generations, Ne_dadi, T, Nu_dadi)
dadi = pd.DataFrame({
    'model': 'dadi',
    'generations': generations,
    'value': dadi_values
})

# Stairway


theta_stair = 4 * Nfinale * mu # 4Nfinale * mu = 4 * (Ninit * ratio_expo) * mu
Ne_stair = theta_stair / (4 * mu * Nu_stair) # N initiale
T = T_stair * 2 * 10000
stair_values = attendu(type_simu, generations, Ne_stair, T, Nu_stair)
stair = pd.DataFrame({
    'model': 'stair',
    'generations': generations,
    'value': stair_values
})



# Regrouper les données
df_all = pd.concat([theo, dadi, stair])

# ------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------- Plot
# ------------------------------------------------------------------------------------
sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.lineplot(data=df_all, x='generations', y='value', hue='model')
plt.xscale('log')
plt.xlabel('Generations (log scale)')
#plt.gca().invert_xaxis()
plt.ylabel('Effective population size (Ne)')
plt.title('Theorical models with infered parameters')
plt.tight_layout()
plt.savefig("../../results/graphs/"  + type_simu  + ".pdf")

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
