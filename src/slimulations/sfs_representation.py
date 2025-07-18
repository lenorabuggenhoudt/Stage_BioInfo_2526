import csv
import argparse
import matplotlib.pyplot as plt
import numpy as np
import warnings
import pandas as pd
import seaborn as sns

# ----------------------
# Argument Parsing
# ----------------------
parser = argparse.ArgumentParser()
parser.add_argument('--replicatenumber', type=int, required=True)
parser.add_argument('--type', type=str, required=True)
args = parser.parse_args()

replicatenumber = args.replicatenumber
typ = args.type

out_dirs = [f"./results/sweep/{typ}", f"./results/neutral/{typ}"]
GR = [1, 10, 50, 100]

# ----------------------
# Data Loading Functions
# ----------------------

# Recup the generation when the simulation was fixed

def gen_fixed(file) :
    with open(file, "r") as f:
        lines = f.readlines()
    gen = 0
    for i, line in enumerate(lines):
        if "FIXED" in line:
            # Find previous generation line
            j = i - 1
            gen_line = lines[j].strip()
            gen = int(gen_line.split(":")[1].strip())
            break
    return gen


#########################
# LOAD ------------------
#########################

data = {
    'GR1': [],
    'GR10': [],
    'GR50': [],
    'GR100': []
}

df_neutral = pd.DataFrame(data)
df_sweep = pd.DataFrame(data)

def load_data(path_prefix, replicatenumber, GR):
    # Temporary structure to accumulate data
    data = {f'GR{i}': [] for i in GR}
    snps = {f'GR{i}': [] for i in GR}
    for r in range(1, replicatenumber + 1):
        for i in GR:
            gr_name = f'GR{i}'
            file = f"{path_prefix}/sfs_{r}_{i}.csv"
            try:
                with open(file, newline='') as csv_file:
                    csv_content = csv.reader(csv_file, delimiter=",")
                    for row in csv_content:
                        # Convert to floats (skip the first element)
                        data_add = [float(x) for x in row[1:]]
                        snp_add = float(row[0])
                        if snp_add < 10 or sum(data_add) > 1.000001 :
                            print("Nbr of snps/sites : ", snp_add)
                            print("Somme data_add : ", sum(data_add))
                            print("R : ", r)
                            print("GR : ", i)
                            print("Path prefix : ", path_prefix)
                        else :
                            snps[gr_name].append(snp_add)
                            data[gr_name].append(data_add)
            except FileNotFoundError:
                print(f"File not found: {file}")
    
    # Create DataFrame from list of lists
    df = pd.DataFrame({k: pd.Series(v) for k, v in data.items()})
    snps = pd.DataFrame({k: pd.Series(v) for k, v in snps.items()})
    return df, snps


df_sweep, snps_sweep = load_data(out_dirs[0], replicatenumber, GR)


# Mean --------------------------------------
def average_sfs_by_gr(df):
    """
    Prend un DataFrame contenant des listes de listes (une par GR)
    et retourne un dictionnaire contenant la moyenne élément par élément
    pour chaque GR.
    """
    averaged_data = {}
    
    for gr_name in df.columns:
        # Récupérer les listes (une par réplica) pour ce GR
        list_of_lists = df[gr_name].dropna().tolist()  # enlever les valeurs manquantes éventuelles

        if list_of_lists:
            # Moyenne élément par élément
            array_data = np.array(list_of_lists)
            averaged_list = np.mean(array_data, axis=0).tolist()
            averaged_data[gr_name] = averaged_list
        else:
            averaged_data[gr_name] = []

    return averaged_data
# ----------------------------------------------------------

averaged_results = average_sfs_by_gr(df_sweep)



# Median -----------------------------------------------------------
def median_sfs_by_gr(df):
    """
    Prend un DataFrame contenant des listes de listes (une par GR)
    et retourne un dictionnaire contenant la médiane élément par élément
    pour chaque GR.
    """
    median_data = {}

    for gr_name in df.columns:
        # Récupérer les listes (une par réplica) pour ce GR
        list_of_lists = df[gr_name].dropna().tolist()

        if list_of_lists:
            # Calcul de la médiane élément par élément
            array_data = np.array(list_of_lists)
            median_list = np.median(array_data, axis=0).tolist()
            median_data[gr_name] = median_list
        else:
            median_data[gr_name] = []

    return median_data

# ------------------------------------------------
median_results = median_sfs_by_gr(df_sweep)

# GR1 neutral ratio -------------------------------------------
# -------------------------------------------------------------


def ratio_gr1_by_gr(df_sweep, df_neutral) :

    average_sfs_sweep = average_sfs_by_gr(df_sweep)
    average_sfs_neutral = average_sfs_by_gr(df_neutral)
    gr1_neutral = average_sfs_neutral['GR1']
    ratios_neutral = {}
    ratios_sweep = {}
    for gr_name, avg in average_sfs_neutral.items():
        if len(avg) == len(gr1_neutral):
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.divide(avg, gr1_neutral)
                ratio = np.nan_to_num(ratio, nan=0.0, posinf=0.0, neginf=0.0)
            ratios_neutral[gr_name] = ratio.tolist()
        else:
            ratios_neutral[gr_name] = []

    for gr_name, avg in average_sfs_sweep.items():
        if len(avg) == len(gr1_neutral):
            with np.errstate(divide='ignore', invalid='ignore'):
                ratio = np.divide(avg, gr1_neutral)
                ratio = np.nan_to_num(ratio, nan=0.0, posinf=0.0, neginf=0.0)
            ratios_sweep[gr_name] = ratio.tolist()
        else:
            ratios_sweep[gr_name] = []

    return ratios_neutral, ratios_sweep



df_sweep, snps_sweep = load_data(out_dirs[0], replicatenumber, GR)
df_neutral, snps_neutral = load_data(out_dirs[1], replicatenumber, GR)

ratios_neutral, ratios_sweep = ratio_gr1_by_gr(df_sweep, df_neutral)


# -------------------------------------------------------------
# -------------------------------------------------------------
# PLOTS 
# -------------------------------------------------------------
# -------------------------------------------------------------

def plot_average_sfs(average_sfs_sweep, average_sfs_neutral, sweep=True, neutral=True):
    """
    Plot average SFS for each GR with same color for sweep and neutral,
    solid line for sweep and dashed line for neutral. No markers.
    """
    gr_names = sorted(average_sfs_sweep.keys())  # to keep order consistent
    num_bins = len(next(iter(average_sfs_sweep.values())))
    half_length = num_bins // 2
    x = list(range(1, half_length + 1))

    
    # Get a color map
    cmap = plt.get_cmap("tab10")
    
    plt.figure(figsize=(12, 6))

    for idx, gr in enumerate(gr_names):
        y_sweep = average_sfs_sweep[gr][:half_length]
        y_neutral = average_sfs_neutral.get(gr, [0]*num_bins)[:half_length]


        color = cmap(idx % 10)  # cycle colors if more than 10 GRs
        if sweep :
            plt.plot(x, y_sweep, label=f'{gr} Sweep', linestyle='--', color=color)
        if neutral :
            plt.plot(x, y_neutral, label=f'{gr} Neutral', linestyle='-', color=color)

    plt.xlabel("Number of individuals with the common mutation")
    plt.ylabel("Percentage of the mutations Having affecting x individuals")
    plt.title("Average Site Frequency Spectrum")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    if not sweep :
        plt.savefig(f"sfs_average_neutral_{typ}.pdf", dpi=300)
    else :
        if not neutral :
            plt.savefig(f"sfs_average_sweep_{typ}.pdf", dpi=300)
        else :
            plt.savefig(f"sfs_average_{typ}.pdf", dpi=300)

# Calculate average SFS for both sweep and neutral
average_sfs_sweep = average_sfs_by_gr(df_sweep)
average_sfs_neutral = average_sfs_by_gr(df_neutral)

plot_average_sfs(average_sfs_sweep, average_sfs_neutral)
plot_average_sfs(average_sfs_sweep, average_sfs_neutral, sweep=False)
plot_average_sfs(average_sfs_sweep, average_sfs_neutral, neutral=False)



# Ratios --------------------------------

def plot_sfs_ratios(ratios_sweep, ratios_neutral, sweep=True, neutral=True):
    """
    Plot only the first half of SFS ratios for each GR with same color.
    Solid line = sweep ratio, dashed line = neutral ratio.
    """
    gr_names = sorted(ratios_sweep.keys())
    full_length = len(next(iter(ratios_sweep.values())))
    half_length = full_length // 2
    x = list(range(1, half_length + 1))

    cmap = plt.get_cmap("tab10")

    plt.figure(figsize=(12, 6))

    for idx, gr in enumerate(gr_names):
        y_sweep = ratios_sweep.get(gr, [0]*full_length)[:half_length]
        y_neutral = ratios_neutral.get(gr, [0]*full_length)[:half_length]
        color = cmap(idx % 10)
        if sweep :
            plt.plot(x, y_sweep, label=f'{gr} Sweep Ratio', linestyle='--', color=color)
        if neutral :
            plt.plot(x, y_neutral, label=f'{gr} Neutral Ratio', linestyle='-', color=color)
    plt.xlabel("Number of individuals with the common mutation")
    plt.ylabel("Ratio (to GR1 Neutral)")
    plt.title("SFS Ratios to GR1 Neutral")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    if not sweep :
        plt.savefig(f"sfs_ratios_neutral_{typ}.pdf", dpi=300)
    else :
        if not neutral :
            plt.savefig(f"sfs_ratios_sweep_{typ}.pdf", dpi=300)
        else :
            plt.savefig(f"sfs_ratios_{typ}.pdf", dpi=300)

ratios_neutral, ratios_sweep = ratio_gr1_by_gr(df_sweep, df_neutral)
plot_sfs_ratios(ratios_sweep, ratios_neutral)
plot_sfs_ratios(ratios_sweep, ratios_neutral, sweep=False)
plot_sfs_ratios(ratios_sweep, ratios_neutral, neutral=False)

#################################
###### BOXPLOT SNPS #############
################################



def prepare_snps_for_boxplot(snps_sweep, snps_neutral):
    data = []
    for gr in snps_sweep.columns:
        # sweep
        vals_sweep = snps_sweep[gr].dropna().tolist()
        data.extend([(gr, 'Sweep', v) for v in vals_sweep])
        # neutral
        vals_neutral = snps_neutral[gr].dropna().tolist()
        data.extend([(gr, 'Neutral', v) for v in vals_neutral])

    df_long = pd.DataFrame(data, columns=['GR', 'Condition', 'SNP'])
    return df_long

# Prépare les données
df_snps_long = prepare_snps_for_boxplot(snps_sweep, snps_neutral)

# Plot
plt.figure(figsize=(10,6))
sns.boxplot(x='GR', y='SNP', hue='Condition', data=df_snps_long, palette='Set2')
plt.title('Distribution of SNPs by GR and Condition')
plt.grid(True)
plt.tight_layout()
plt.savefig(f'boxplots_snps_{typ}.pdf', dpi=300)