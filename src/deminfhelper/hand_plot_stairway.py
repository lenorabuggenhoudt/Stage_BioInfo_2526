import numpy as np
import matplotlib.pyplot as plt
import os
import argparse
import sys
import time
# -------------------------------------------
# ------------------ Arguments --------------
# -------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument('--config_file', type=str, required=True)
parser.add_argument('--output_dir', type=str, required=True)
parser.add_argument('--T_theo', type=float, required=True)
parser.add_argument('--xlog', type=int, required=True)
parser.add_argument('--ylog', type=int, required=True)
parser.add_argument('--GR', type=int, required=True)
args = parser.parse_args()

config_file = args.config_file # Path to the Stairway output file of DFH
output_dir = args.output_dir
T_theo = args.T_theo
x_log = bool(args.xlog)
y_log = bool(args.ylog)
# -------------------------------------------
import yaml
with open(config_file, 'r') as f:
    data = yaml.full_load(f)
stairway_file = data.get("out_dir_stairwayplot2") + "cerevisiae/cerevisiae.final.summary"
if 'expansion' in data.get('vcf') :
    simu_type = 'expansion_'
else :
    simu_type = 'reduction_'
if 'neutral' in data.get('vcf') :
    simu_type = simu_type + 'neutral'
else :
    simu_type = simu_type + 'sweep'

output_file = output_dir + data.get('name_pop') 

# --------------------------------------------

def plot_stairwayplot2(popid, summary_file, out_dir, T_theo, xlog=True, ylog=True):
    """
    Generate a stairway plot from summary data.

    This function reads summary data from a Stairway Plot analysis, extracts the
    effective population size (Ne) estimates over time, and generates a plot.

    Parameters:
        popid (str): Identifier for the population.
        summary_file (str): Path to the Stairway Plot summary data file.
        out_dir (str): Directory where the stairway plot image will be saved.

    Note:
        - The plot will be saved as "{popid}_stw_plot.png" in the specified 'out_dir'.

    Returns:
        None: The function saves the stairway plot as an image file but does not return any values.
    """
    Ne_med=[]
    Ne1=[]
    Ne3=[]
    T=[]
    with open(summary_file) as input_file:
        line = input_file.readline()
        line = input_file.readline()
        while line != '':
            #if float(line.split('\t')[6]) not in Ne or float(line.split('\t')[5]) not in T:
            Ne_med.append(float(line.split('\t')[6]))
            Ne1.append(float(line.split('\t')[7]))
            Ne3.append(float(line.split('\t')[8]))
            T.append(float(line.split('\t')[5]))
            line = input_file.readline()
    plt.plot(T,Ne_med,color="red", marker="x")
    plt.plot(T,Ne1,color="grey")
    plt.plot(T,Ne3,color="grey")

    #theorical values
    plt.vlines(x=[T_theo], ymin=0, ymax=max(Ne_med), colors=['g'], linestyles=['--'])


    plt.title(fr"[StairwayPlot2] {popid}")
    logx = ''
    logy =''
    if xlog :
        plt.xscale("log")
        logx ="_logx"
    if ylog:
        plt.yscale("log")
        logy ="_logy"
    plt.xlabel('Time (in years)')
    plt.ylabel('Ne')
    plt.savefig(output_file + logx + logy + "_" + simu_type + "_" +  str(data.get('GR')) + "_stw_plot_" + str(time.time()) + ".pdf")
    plt.savefig(data.get('out_dir_stairwayplot2') + data.get('name_pop') + "_" + logx + logy + "_" + simu_type + "_" + str(data.get('GR')) + "_stw_plot_" + ".pdf")
    plt.close()
print("XLOG : ", x_log)
print("YLOG : ", y_log)
plot_stairwayplot2('1', stairway_file, output_file, T_theo, xlog=x_log, ylog=y_log )
if not (x_log and y_log) :
    plot_stairwayplot2('1', stairway_file, output_file, T_theo, xlog=True, ylog=True)
