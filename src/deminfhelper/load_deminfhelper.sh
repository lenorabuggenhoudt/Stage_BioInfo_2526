#!/bin/bash
                                                                                     
## Lab : LISN ~ UMR9015 ~ BIOINFO team

#SBATCH --job-name=deminfhelper_launch
#SBATCH --array=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --account=inferscerevol

# Necessary libraries 
# added
module load java-jdk/8.0.112
module load singularity
module load seqtk

module load bcftools/1.16 
module load dadi/2.2.0 

module load mysql-client/8.0.23 
module load psmc/0.6.5 
module load smcpp
module load vcftools/0.1.16

module load python/3.12 #python/3.10 # test to try and debug smc++
pip install numpy==1.23.5

# ---------------- recup args --------------------------
for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   export "$KEY"="$VALUE"
done
#out_dir=$OUT_DIR
#echo "OUT_DIR = $OUT_DIR"
#vcf_dir=$VCF_DIR
#echo "VCF_DIR = $VCF_DIR"

# ---------------- Optional arguments (CAN be given by the user but have default values) ----------------

if [ -v GEN_TIME ]; then
	echo "GEN_TIME = $GEN_TIME"
	gen_time=$GEN_TIME # Population size 
else 
	gen_time=0.0027 # (1 per day for a year) 
fi

if [ -v MUT_RATE ]; then
	echo "MUT_RATE = $MUT_RATE"
	mut_rate=$MUT_RATE # Population size 
else 
	mut_rate=1e-7 # (1 per day for a year) 
fi

if [ -v P0 ]; then
	echo "P0 = $P0"
	p0=$P0 # Population size 
else 
	p0="1, 1" 
fi

if [ -v LBOUND ]; then
	echo "LBOUND = $LBOUND"
	lbound=$LBOUND # Population size 
else 
	lbound="0.1, 0.1" 
fi

if [ -v UBOUND ]; then
	echo "UBOUND = $UBOUND"
	ubound=$UBOUND # Population size 
else 
	ubound="100, 100" 
fi

# creation fichier config
#config_file="${outdir}config_file.yml"
#python3 config_creation.py --out_dir $outdir --vcf_dir $vcfdir --generation_time $gen_time --mutation_rate $mut_rate --p0 $p0 --lower_bound $lbound --upper_bound $ubound > $config_file





# ----------------------------------------------------


config_file=$1
echo "$config_file"


# ----------------------- DFH ------------------------
# calls to deminfhelper
#echo "Appel SFS"
#python3 ../deminfhelper/deminfhelper.py --config_file $config_file --sfs > Sim.out
#echo "Appel Dadi"
#dadi-cli Model --names 
#python3 ../deminfhelper/deminfhelper.py --config_file  $config_file --dadi

#echo "Appel MSMC2"
#python3 ../deminfhelper/deminfhelper.py --config_file $config_file --msmc2 --cpus 10
#python3 ../deminfhelper/deminfhelper.py --config_file $config_file --plot_msmc2 --cpus 10
# test mail TF
#python3 ../deminfhelper/deminfhelper.py --config_file $config_file --msmc2 --plot_msmc2 --cpus 10

#echo "Appel PSMC"
#python3 ../deminfhelper/deminfhelper.py --config_file $config_file --psmc --cpus 10
#python3 ../deminfhelper/deminfhelper.py --config_file $config_file --plot_psmc
# Swp2
#echo "Appel Stairway plot"
#python3 ../deminfhelper/deminfhelper.py --config_file $config_file --stairwayplot2
#python3 ../deminfhelper/hand_plot_stairway.py --config_file $config_file --output_dir '../../results/stairway_plots/' --xlog 0 --ylog 0 --T_theo 5.48

# SMC++: problème de ##contig=<ID=1,length=2200000>, longueur des contigs qui ne match pas dans le header avec les noms des contigs "chr1,chr2..."
#echo "SMC ++"
#echo "Avec un module load de Python 3.10"
python3 ../deminfhelper/deminfhelper.py --config_file $config_file --smcpp
