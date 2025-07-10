import csv
import argparse
import matplotlib.pyplot as plt
import time

parser = argparse.ArgumentParser()
parser.add_argument('--replicatenumber', type=int, required=True)
parser.add_argument('--generationmitoses', type=int, required=True)
parser.add_argument('--out_dir', type=str, required=True)
args = parser.parse_args()

replicatenumber = args.replicatenumber
generationmitoses = args.generationmitoses
out_dir = args.out_dir

file = f"./{out_dir}/sfs_{replicatenumber}_{generationmitoses}.csv"

with open(file, newline='') as csv_file:
    csv_content = csv.reader(csv_file, delimiter=",")
    data = []
    for ligne in csv_content:
        # ligne is a list of strings from one row; assuming your file has a single row of comma-separated floats
        data = [float(x.strip()) for x in ligne]


# -----------------------------------------------------------------------------
# --------------------- HISTOGRAM ---------------------------------------------
# -----------------------------------------------------------------------------


x = list(range(len(data)))

fig = plt.figure(figsize=(12, 6))
plt.bar(x, data, width=1.0, color='steelblue')

plt.xlabel('Number of site in commun')
plt.ylabel('Frequence (echelle * 10000)')
plt.title('Histogram Allele Frequency Spectrum')

ax = fig.gca()
plt.setp(ax.get_xticklabels(), visible=False)       # hide all x tick labels
plt.setp(ax.get_xticklabels()[::5], visible=True)   # show every 5th x tick label

plt.tight_layout()
plt.savefig(f"./{out_dir}/histo_sfs_{replicatenumber}_{generationmitoses}_{time.time()}.pdf")


