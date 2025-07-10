#!/bin/bash
                                                                                     
## Lab : LISN ~ UMR9015 ~ BIOINFO team

#SBATCH --job-name=softwares_launch
#SBATCH --array=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --account=inferscerevol

# Necessary libraries 
# added


# Necessary libraries 
# added
module load java-jdk/8.0.112
module load singularity
module load seqtk

module load bcftools/1.16 
module load dadi/2.2.0 

module load mysql-client/8.0.23 
module load psmc/0.6.5 
module load smcpp/1.15.4 
module load vcftools/0.1.16

module load python/3.12 #python/3.10 # test to try and debug smc++

# ---------------- recup args --------------------------
for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   export "$KEY"="$VALUE"
done

# Mandatory arguments

echo "TYPE = $TYPE" # 'reduction' or 'expansion'
echo "SWEEP = $SWEEP" # 'true' or 'false'
if [ $SWEEP == "false" ]; then  
	sweep='neutral'
else 
	sweep='sweep'
fi
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

if [ -v T_THEO ]; then
	echo "T_Theo = $T_THEO"
	T_Theo=$T_THEO # Theorical T
else 
	T_Theo=0.1
fi

if [ -v GR ]; then
	echo "GR = $GR"
	GR=$GR # Theorical T
else 
	GR=1
fi

if [ -v REPET ]; then
	echo "REPET = $REPET"
	repetition_number=$REPET # Theorical T
else 
	repetition_number=1
fi



# ----------------------- Loop ------------------------
for generationmitoses in 10; do

	config_file="${TYPE}_${sweep}_${generationmitoses}_config.yml"
	echo "Config file : $config_file "
	# calls to deminfhelper
	#echo "Appel SFS"
	#python3 ../deminfhelper/deminfhelper.py --config_file $config_file --sfs > Sim.out

	# Stairway
	echo "Stairway"
	python3 ../deminfhelper/deminfhelper.py --config_file $config_file --stairwayplot2
	#echo "Stairway plot "
	#python3 ../deminfhelper/hand_plot_stairway.py --config_file $config_file --output_dir '../../results/stairway_plots/' --xlog 0 --ylog 0 --T_theo $T_Theo --GR $generationmitoses

	# Dadi
	#echo "Dadi "
	#python3 ../dadi_louis/script_dadi_louis.py --config_file $config_file --GR $generationmitoses --repetiton_number $repetition_number


	done
