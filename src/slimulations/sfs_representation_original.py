import csv
import argparse
import matplotlib.pyplot as plt
import numpy as np
import warnings

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
    print("Gen renvoyé : ", gen)
    return gen





def load_data(path_prefix, replicatenumber, GR):
    data_all = []
    snps_all = []
    gen_all = []
    for r in range(1, replicatenumber + 1):
        data_r = []
        snps_r = []
        gen_r = []
        for i, gr in enumerate(GR):
            file = f"{path_prefix}/sfs_{r}_{gr}.csv"
            runinfos = f"{path_prefix}/runInfo_{r}_{gr}.txt"
            gen = gen_fixed(runinfos)
            gen_r.append(gen) 
            with open(file, newline='') as csv_file:
                csv_content = csv.reader(csv_file, delimiter=",")
                for ligne in csv_content:
                    data_add = [float(x) for x in ligne[1:]]
                    snps_add = float(ligne[0])
                    data_r.append(data_add)
                    snps_r.append(snps_add)
        data_all.append(data_r)
        snps_all.append(snps_r)
        gen_all.append(gen_r)
    return data_all, snps_all, gen_all

data_sweep, snps_sweep, gen_sweep = load_data(out_dirs[0], replicatenumber, GR)
data_neutral, snps_neutral, gen_neutral = load_data(out_dirs[1], replicatenumber, GR)
print("Data : sweep : ", data_sweep)
# ----------------------
# Mean Frequency
# ----------------------
def get_mean_freq(data, snps, GR_index):
    print( "data : ", data)
    n_classes = len(data[0][GR_index])
    rep = len(data)
    mean_data = [0.0] * n_classes
    for i in range(rep):
        for j in range(n_classes):
            mean_data[j] += data[i][GR_index][j] 
    mean_data = [x / rep for x in mean_data]
    
    # check egal to 1
    summ = 0
    for el in mean_data :
        summ += el
    print("Sum :", summ)
    print('GR : ', GR_index)
    print("MD : ", mean_data)
    return mean_data

# -----------------------
# Median Frequency
# -----------------------
import statistics

def get_median_freq(data, snps, GR_index):
    n_classes = len(data[0][GR_index])
    class_values = [[] for _ in range(n_classes)]

    for i in range(len(data)):
        if snps[i][GR_index] == 0:
            continue
        for j in range(n_classes):
            # Normalize the frequency before appending
            value = data[i][GR_index][j] / snps[i][GR_index]
            class_values[j].append(value)

    # Compute median for each class
    median_data = [
        statistics.median(values) if values else 0.0
        for values in class_values
    ]
    return median_data



# --------------------
# Get the snps
# --------------------
def get_snps(snps, GR) :
    res = []
    for i in range(len(snps)): # rep
        res.append(snps[i][GR])
    return res

def get_genfixed(gen, GR) :
    res = []
    for i in range(len(gen)): # rep
        res.append(gen[i][GR])
    return res
# ----------------------
#Division Function
# ----------------------
def safe_divide(numerator, denominator):
    numerator = np.array(numerator)
    denominator = np.array(denominator)
    with np.errstate(divide='ignore', invalid='ignore'):
        result = np.true_divide(numerator, denominator)
        result[~np.isfinite(result)] = 0.0  # replace inf, -inf, nan with 0
    return result

# ----------------------
# Compute Frequencies and Ratios
# ----------------------
print("Size : ", len(data_neutral[0][0]))
size_kept = int(len(data_neutral[0][0]) / 2)
with np.errstate(invalid='ignore'):
    GR1_neutral = get_mean_freq(data_neutral, snps_neutral, 0)[0:size_kept]
    GR1_sweep = get_mean_freq(data_sweep, snps_sweep, 0)[0:size_kept]

    GR1_neutral_med = get_median_freq(data_neutral, snps_neutral, 0)[0:size_kept]
    GR1_sweep_med = get_median_freq(data_sweep, snps_neutral, 0)[0:size_kept]

    GR1_neutral_ratio = safe_divide(GR1_neutral_med, GR1_neutral_med)
    GR1_sweep_ratio = safe_divide(GR1_sweep_med, GR1_neutral_med)
    neutral_gr1_snps = get_snps(snps_neutral, 0)
    sweep_gr1_snps = get_snps(snps_sweep, 0)
    neutral_gr1_genfixed = get_genfixed(gen_neutral, 0)
    sweep_gr1_genfixed = get_genfixed(gen_sweep, 0)

# -----
    GR10_neutral = get_mean_freq(data_neutral, snps_neutral, 1)[0:size_kept]
    GR10_sweep = get_mean_freq(data_sweep, snps_sweep, 1)[0:size_kept]

    GR10_neutral_med = get_median_freq(data_neutral, snps_neutral, 1)[0:size_kept]
    GR10_sweep_med = get_median_freq(data_sweep, snps_neutral, 1)[0:size_kept]

    GR10_neutral_ratio = safe_divide(GR10_neutral_med, GR1_neutral_med)
    GR10_sweep_ratio = safe_divide(GR10_sweep_med, GR1_neutral_med)
    neutral_gr10_snps = get_snps(snps_neutral, 1)
    sweep_gr10_snps = get_snps(snps_sweep, 1)
    neutral_gr10_genfixed = get_genfixed(gen_neutral, 1)
    sweep_gr10_genfixed = get_genfixed(gen_sweep, 1)
# -----
    GR50_neutral = get_mean_freq(data_neutral, snps_neutral, 2)[0:size_kept]
    GR50_sweep = get_mean_freq(data_sweep, snps_sweep, 2)[0:size_kept]

    GR50_neutral_med = get_median_freq(data_neutral, snps_neutral, 2)[0:size_kept]
    GR50_sweep_med = get_median_freq(data_sweep, snps_neutral, 2)[0:size_kept]

    GR50_neutral_ratio = safe_divide(GR50_neutral_med, GR1_neutral_med)
    GR50_sweep_ratio = safe_divide(GR50_sweep_med, GR1_neutral_med)
    neutral_gr50_snps = get_snps(snps_neutral, 2)
    sweep_gr50_snps = get_snps(snps_sweep, 2)
    neutral_gr50_genfixed = get_genfixed(gen_neutral, 2)
    sweep_gr50_genfixed = get_genfixed(gen_sweep, 2)
# ----
    GR100_neutral = get_mean_freq(data_neutral, snps_neutral, 3)[0:size_kept]
    GR100_sweep = get_mean_freq(data_sweep, snps_sweep, 3)[0:size_kept]

    GR100_neutral_med = get_median_freq(data_neutral, snps_neutral, 3)[0:size_kept]
    GR100_sweep_med = get_median_freq(data_sweep, snps_neutral, 3)[0:size_kept]

    GR100_neutral_ratio = safe_divide(GR100_neutral_med, GR1_neutral_med)
    GR100_sweep_ratio = safe_divide(GR100_sweep_med, GR1_neutral_med)
    neutral_gr100_snps = get_snps(snps_neutral, 3)
    sweep_gr100_snps = get_snps(snps_sweep, 3)
    neutral_gr100_genfixed = get_genfixed(gen_neutral, 3)
    sweep_gr100_genfixed = get_genfixed(gen_sweep, 3)
    


# ----------------------
# Plotting
# ----------------------
x = list(range(1, len(GR1_neutral_ratio) + 1))

# --- Combined plot ---
plt.figure(figsize=(10, 6))
plt.plot(x, GR1_neutral_ratio, label="GR=1 (neutral)", color="blue")
plt.plot(x, GR10_neutral_ratio, label="GR=10 (neutral)", color="green")
plt.plot(x, GR50_neutral_ratio, label="GR=50 (neutral)", color="orange")
plt.plot(x, GR100_neutral_ratio, label="GR=100 (neutral)", color="purple")
plt.plot(x, GR1_sweep_ratio, label="GR=1 (sweep)", linestyle="-.", color="blue")
plt.plot(x, GR10_sweep_ratio, label="GR=10 (sweep)", linestyle="-.", color="green")
plt.plot(x, GR50_sweep_ratio, label="GR=50 (sweep)", linestyle="-.", color="orange")
plt.plot(x, GR100_sweep_ratio, label="GR=100 (sweep)", linestyle="-.", color="purple")
plt.xlabel("Number of sites in common")
plt.ylabel("Relative SFS (normalized by GR=1 neutral)")
plt.title("SFS Ratio per Generation Rate (Neutral vs Sweep)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"sfs_comparison_{typ}.pdf", dpi=300)

# --- Sweep-only plot ---
plt.figure(figsize=(10, 6))
plt.plot(x, GR1_neutral_ratio, label="GR=1 (neutral)", color="blue")
plt.plot(x, GR1_sweep_ratio, label="GR=1 (sweep)", linestyle="-.", color="blue")
plt.plot(x, GR10_sweep_ratio, label="GR=10 (sweep)", linestyle="-.", color="green")
plt.plot(x, GR50_sweep_ratio, label="GR=50 (sweep)", linestyle="-.", color="orange")
plt.plot(x, GR100_sweep_ratio, label="GR=100 (sweep)", linestyle="-.", color="purple")
plt.xlabel("Number of sites in common")
plt.ylabel("Relative SFS (normalized by GR=1 neutral)")
plt.title("SFS Ratio per Generation Rate (Sweep)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"sfs_comparison_{typ}_sweep.pdf", dpi=300)

# --- Neutral-only plot ---
plt.figure(figsize=(10, 6))
plt.plot(x, GR1_neutral_ratio, label="GR=1 (neutral)", color="blue")
plt.plot(x, GR10_neutral_ratio, label="GR=10 (neutral)", color="green")
plt.plot(x, GR50_neutral_ratio, label="GR=50 (neutral)", color="orange")
plt.plot(x, GR100_neutral_ratio, label="GR=100 (neutral)", color="purple")
plt.xlabel("Number of sites in common")
plt.ylabel("Relative SFS (normalized by GR=1 neutral)")
plt.title("SFS Ratio per Generation Rate (Neutral)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"sfs_comparison_{typ}_neutral.pdf", dpi=300)


# --- SNPS  ---
groups = ['GR1', 'GR10', 'GR50', 'GR100']

neutral_data = [neutral_gr1_snps, neutral_gr10_snps, neutral_gr50_snps, neutral_gr100_snps]
sweep_data = [sweep_gr1_snps, sweep_gr10_snps, sweep_gr50_snps, sweep_gr100_snps]

# Position des boxplots sur l'axe x
positions_neutral = [i*2 for i in range(len(groups))]      # 0, 2, 4, 6
positions_sweep = [i*2 + 1 for i in range(len(groups))]    # 1, 3, 5, 7

plt.figure(figsize=(10,6))

# Boxplots neutral (gauche)
plt.boxplot(neutral_data, positions=positions_neutral, widths=0.6, patch_artist=True,
            boxprops=dict(facecolor="lightblue"), medianprops=dict(color="blue"), labels=None)

# Boxplots sweep (droite)
plt.boxplot(sweep_data, positions=positions_sweep, widths=0.6, patch_artist=True,
            boxprops=dict(facecolor="lightcoral"), medianprops=dict(color="red"), labels=None)

# X axis ticks au centre des groupes
xticks_positions = [(n + s)/2 for n, s in zip(positions_neutral, positions_sweep)]
plt.xticks(xticks_positions, groups)

plt.xlabel("Groupes GR")
plt.ylabel("Nombre de SNPs")
plt.title("Comparaison du nombre des SNPs par groupe et condition")
plt.plot([], c="lightblue", label="Neutral")
plt.plot([], c="lightcoral", label="Sweep")
plt.legend()

plt.savefig(f"snps_{typ}.pdf", dpi=300)

# ------------------------
# Gen time fixed
# -----------------------

groups = ['GR1', 'GR10', 'GR50', 'GR100']
sweep_data = [sweep_gr1_genfixed, sweep_gr10_genfixed, sweep_gr50_genfixed, sweep_gr100_genfixed]

plt.figure(figsize=(10,6))

plt.boxplot(sweep_data)

xticks_positions = [r for r in range(1, len(sweep_data) + 1)]
plt.xticks(xticks_positions, groups)

plt.xlabel("Groupes GR")
plt.ylabel("Fixation generation")
plt.title("Comparaison du temps de fixation de la mutation par groupe et condition")
plt.plot([], c="lightblue", label="Neutral")
plt.plot([], c="lightcoral", label="Sweep")
plt.legend()

plt.savefig(f"generation_fixed_{typ}.pdf", dpi=300)

# ---------------------
# Plot w/o ratio
# ---------------------

plt.figure(figsize=(10, 6))
plt.plot(x, GR1_neutral, label="GR=1 (neutral)", color="blue")
plt.plot(x, GR10_neutral, label="GR=10 (neutral)", color="green")
plt.plot(x, GR50_neutral, label="GR=50 (neutral)", color="orange")
plt.plot(x, GR100_neutral, label="GR=100 (neutral)", color="purple")
plt.plot(x, GR1_sweep, label="GR=1 (sweep)", linestyle="-.", color="blue")
plt.plot(x, GR10_sweep, label="GR=10 (sweep)", linestyle="-.", color="green")
plt.plot(x, GR50_sweep, label="GR=50 (sweep)", linestyle="-.", color="orange")
plt.plot(x, GR100_sweep, label="GR=100 (sweep)", linestyle="-.", color="purple")
plt.xlabel("Number of individuals concerned")
plt.ylabel("Number of mutations")
plt.title("Mean SFS per Generation Rate (Neutral vs Sweep)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"sfs_{typ}.pdf", dpi=300)

# ---------------------
# Plot w/o ratio et median
# ---------------------

plt.figure(figsize=(10, 6))
plt.plot(x, GR1_neutral_med, label="GR=1 (neutral)", color="blue")
plt.plot(x, GR10_neutral_med, label="GR=10 (neutral)", color="green")
plt.plot(x, GR50_neutral_med, label="GR=50 (neutral)", color="orange")
plt.plot(x, GR100_neutral_med, label="GR=100 (neutral)", color="purple")

plt.plot(x, GR1_sweep_med, label="GR=1 (sweep)", linestyle="-.", color="blue")
plt.plot(x, GR10_sweep_med, label="GR=10 (sweep)", linestyle="-.", color="green")
plt.plot(x, GR50_sweep_med, label="GR=50 (sweep)", linestyle="-.", color="orange")
plt.plot(x, GR100_sweep, label="GR=100 (sweep)", linestyle="-.", color="purple")

plt.xlabel("Number of individuals concerned")
plt.ylabel("Number of mutations")
plt.title("SFS median per Generation Rate (Neutral vs Sweep)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(f"sfs_{typ}_median.pdf", dpi=300)

