import argparse
import matplotlib.pyplot as plt
import subprocess
import pyslim
import tskit
import os
import random

parser = argparse.ArgumentParser()
parser.add_argument('--replicatenumber', type=int, required=True)
parser.add_argument('--recombinationrate', type=float, required=True)
parser.add_argument('--generationmitoses', type=int, required=True)
parser.add_argument('--mutationrate', type=float, required=True)
parser.add_argument('--sweep_coeff', type=float, required=True)
parser.add_argument('--dom_coeff', type=float, required=True)
parser.add_argument('--reduc_ratio', type=float, required=True)
parser.add_argument('--exp_ratio', type=float, required=True)
parser.add_argument('--out_dir', type=str, required=True)
parser.add_argument('--slim_file', type=str, required=True)
parser.add_argument('--trees_file', type=str, required=True)
parser.add_argument('--popsize', type=str, required=True)
args = parser.parse_args()



rep = args.replicatenumber
recombination_rate = args.recombinationrate
generationmitoses = args.generationmitoses
mutation_rate = args.mutationrate
sweep_coeff = args.sweep_coeff
dom_coeff = args.dom_coeff
reduc_ratio = args.reduc_ratio
exp_ratio = args.exp_ratio
out_dir = args.out_dir
slim_file = args.slim_file
trees_file = args.trees_file
popsize = args.popsize


mutation_site = 500000

trees_file = f"{out_dir}/ts_{rep}_{generationmitoses}.trees"
tmpSLiM=f"./{out_dir}/runInfo_{rep}_{generationmitoses}.txt" #  Temporary file, will contain informations about the run (seed, time, etc)
print("Tree file ", trees_file)
# Function to run SLiM with a new seed
def run_slim_once(seed):
    result = subprocess.run(
        ["slim", "-d", f"popsize={popsize} ", 
                "-d", f"recombinationrate={recombination_rate} ", 
                "-d", f"replicatenumber={rep}", 
                "-d", f"generationmitoses={generationmitoses}",
                "-d", f"mutation_rate={mutation_rate} ",
                "-d", f"sweep_coeff={sweep_coeff} ",
                "-d", f"dominance_coeff={dom_coeff} ",
                "-d", f"reductionratio={reduc_ratio} ",
                "-d", f"expansionratio={exp_ratio} ",
                "-d", f"seed={seed} ",
                slim_file
        ],  
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    return result.returncode == 0

# Function to check if sweep fixed
def check_sweep_status(ts, sweep_site):
    sweep_muts = [m for m in ts.mutations() if m.site == sweep_site and ts.site(m.site).position == sweep_site]
    
    if len(sweep_muts) == 0:
        return "lost"
    
    # Check if it's fixed (present in all individuals)
    derived_node = sweep_muts[0].node
    is_fixed = all(derived_node in tree.nodes() for tree in ts.trees())
    return "fixed" if is_fixed else "not fixed"

# Retry loop
max_attempts = 100
for attempt in range(max_attempts):
    print(f"\nAttempt {attempt + 1}")

    # Run SLiM
    seed = random.randint(1, 2**62 - 1)
    success = run_slim_once(seed)

    # Load tree sequence
    if not os.path.exists(trees_file):
        print("Tree file missing after SLiM run.")
        continue

    ts = tskit.load(trees_file)

    # Check sweep status
    status = check_sweep_status(ts, mutation_site)
    print(f"Sweep status: {status}")

    if status != "lost":
        print("Sweep not lost, keep this replicate.")
        break
    else:
        print("Sweep lost, discarding and retrying...")
        os.remove(trees_file)

else:
    print("Max attempts reached without fixation.")