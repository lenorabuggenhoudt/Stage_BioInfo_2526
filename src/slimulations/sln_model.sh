#!/bin/bash
                                                                                     
## Lab : LISN ~ UMR9015 ~ BIOINFO team

#SBATCH --job-name=slimulations
#SBATCH --array=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --account=inferscerevol

# -------------------------- Goal --------------------------

## This script will launch SLiM simulations that will generate some tree sequences, then 
## use Python to calculate some sumstats for the given model,
## and finally create a vcf file for each simulation.

# -----------------------------------------------------------------------

# -------------------------- Environment setup --------------------------

# -----------------------------------------------------------------------

# Comment the conda activate (we do not use conda)
#source /shared/ifbstor1/software/miniconda/bin/activate SLiMLouis
# Beware that the module load order seems to have an effect 
module load vcftools
module load htslib   # to be able to call bgzip and tabix
module load slim/4.0.1
module load tskit/0.5.3
module load msprime/1.2.0
module load scikit-allel
module load pyslim

module load python/3.12

# -----------------------------------------------------------------------------------------

# -------------------------- Parameters for the simulations -------------------------------

# -----------------------------------------------------------------------------------------

# Necessary arguments :
# -  TYPE = "reduction" or "expansion", represents the type of simulation wanted
# -  SWEEP = "false" or "true", represents the presence or absence of selection during the simulation ((if mispelled considered as true))

# Optional arguments 
# -  POPSIZE = int, represents the initial size of the population we will simulate (by default is 10000)
# -  NB_SAMPLES = int, the number of samples (individuals) we will keep in the VCF file (by default is 100)  
#				/!\ Beware that NB_SAMPLES needs to be lesser than the FINAL size of the population (after the simulation)

# -  MUT_RATE = float, the mutation rate we will apply in the simulation (by default is 1e-7)
# -  KEEP_VCF = "false" or "true", indicates if we want to conserve the original VCF created (without any of the handmade modifications, that would in that case be saved in another file (.modified))
#               (by default is considered false)
# -  B_RATIO = int (>=0), indicates the reduction ratio 
# -  E_RATIO = int (>=0), indicates the expansion ratio

# Example of working calls :
# sbatch sln_model.sh TYPE="reduction" SWEEP="false"
# sbatch sln_model.sh TYPE="reduction" SWEEP="false" POPSIZE=1000
# sbatch sln_model.sh TYPE="reduction" SWEEP="false" POPSIZE=100 NB_SAMPLES=5
# sbatch sln_model.sh TYPE="reduction" SWEEP="false" MUT_RATE=1e-9 KEEP_VCF="true" POPSIZE=1000


# ---------------- Necessary arguments (given by the user in the command line) ----------------

# We gather the arguments given in the command line
for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   export "$KEY"="$VALUE"
done
echo "TYPE = $TYPE"
echo "SWEEP = $SWEEP"

# ---------------- Optional arguments (CAN be given by the user but have default values) ----------------

if [ -v POPSIZE ]; then
	echo "POPSIZE = $POPSIZE"
	popsize=$POPSIZE # Population size 
else 
	popsize=10000 # Population size 
fi

if [ -v NB_SAMPLES ]; then
	echo "NB_SAMPLES = $NB_SAMPLES"
	nb_samples=$NB_SAMPLES # Number of samples for the python part
else 
	nb_samples=100 # Number of samples for the python part
fi

if [ -v MUT_RATE ]; then
	echo "MUT_RATE = $MUT_RATE"
	mutation_rate=$MUT_RATE # Mutation rate  
else 
	mutation_rate=1e-7 # Mutation rate 
fi

if [ -v KEEP_VCF ]; then
	echo "KEEP_VCF = $KEEP_VCF"
	keep_vcf=$KEEP_VCF # indicates if we want to keep the original VCF file (our modifications will be saved on another VCF file (.modified.vcf))  
else 
	keep_vcf="false" 
fi

if [ -v B_RATIO ]; then
	echo "B_RATIO = $B_RATIO"
	reductionratio=$B_RATIO # indicates the reduction ratio (1/reductionratio * taille population initiale = taille population finale)  
else 
	reductionratio=10 # By default this ratio is 10 
fi

if [ -v E_RATIO ]; then
	echo "E_RATIO = $E_RATIO"
	expansionratio=$E_RATIO # indicates the expansion ratio (expansionratio = taille population finale / taille population initiale)    
else 
	expansionratio=15 # By default this ratio is 15
fi


if [ -v REP ]; then
	echo "REP = $REP"
	replicatenumber=$REP # indicates the expansion ratio (expansionratio = taille population finale / taille population initiale)    
else 
	replicatenumber=1 # By default this ratio is 15
fi


if [ -v ALL ]; then
	echo "ALL = $ALL"
	all=$ALL # indicates the expansion ratio (expansionratio = taille population finale / taille population initiale)    
else 
	all='true' # By default this ratio is 15
fi

if [ -v ADD ]; then
	echo "ADD = $add"
	add=$ADD # indicates the expansion ratio (expansionratio = taille population finale / taille population initiale)    
else 
	add=$1 # By default this ratio is 15
fi
# ---------------- Hard coded arguments ----------------

windows=50 # Number of windows for sumstats computation on the first chromosome
recombinationrate=1e-7 # Recombination rate (not scaled with GR!) (rho)
sweep_coeff=0.05 # Selection coefficient
dominance_coeff=0.5 # Dominance coefficient 
#replicatenumber=${SLURM_ARRAY_TASK_ID} # Replicate number


# ----------------------------------------------------

# ------------------ Simulations ---------------------

# ----------------------------------------------------
# Run the SLiM simulation : generate a .trees file (= a tree sequence)
# We select the right slim file and the right output directory
if [[ $SWEEP == "false" ]]; then
	echo "neutral"
	
	if [ $TYPE == "reduction" ]; then
		echo "reduction"
		slim_file="reduction_neutral.slim"
		out_dir="results/neutral/reduction"
	else
		if [ $TYPE == "expansion" ]; then
			echo "expansion"
			slim_file="expansion_neutral.slim"
			out_dir="results/neutral/expansion"
		else
			echo "constant"
			slim_file="cst_neutral.slim"
			out_dir="results/neutral/constant"
		fi
		
	fi
else
	echo "sweep"
	if [ $TYPE == "reduction" ]; then
		echo "reduction"
		slim_file="reduction_sweep.slim"
		out_dir="results/sweep/reduction"
	else
		if [ $TYPE == "expansion" ]; then
			echo "expansion"
			slim_file="expansion_sweep.slim"
			out_dir="results/sweep/expansion"
		else
			echo "constant"
			slim_file="cst_sweep.slim"
			out_dir="results/sweep/constant"
		fi
	fi
fi # to end the if statement

if [ "$all" == "true" ]; then
    start=$replicatenumber
    fin=$(( replicatenumber + add ))
    echo "fin : $fin"
else
    start=$replicatenumber
    fin=$replicatenumber
    echo "fin : $fin"
fi

for repet in $(seq $start $fin )
do
# Test different value for generationmitoses (= 1/alpha) 
	echo "r : $repet"
	for generationmitoses in 1 10 50 100; do

		tmpSLiM="./${out_dir}/runInfo_${repet}_${generationmitoses}.txt" #  Temporary file, will contain informations about the run (seed, time, etc)

		slim -d popsize=$popsize -d recombinationrate=$recombinationrate -d replicatenumber=$repet -d generationmitoses=$generationmitoses -d mutation_rate=$mutation_rate -d sweep_coeff=$sweep_coeff -d dominance_coeff=$dominance_coeff -d reductionratio=$reductionratio -d expansionratio=$expansionratio $slim_file > $tmpSLiM

		trees_file="${out_dir}/ts_${repet}_${generationmitoses}.trees"
		
		# Run the Python script : compute SFS, do the demography parameters inference and save them in a file
		#python3 tree_sumstats.py --replicatenumber $repet --nb_samples $nb_samples --popsize $popsize --mutation_rate $mutation_rate --recombinationrate $recombinationrate --windows $windows --generationmitoses $generationmitoses --out_dir $out_dir

		# Transform the trees file into a vcf file
		python3 tree_to_vcf.py --replicatenumber $repet --generationmitoses $generationmitoses --recombinationrate $recombinationrate --mutation_rate $mutation_rate --popsize $popsize --nb_samples $nb_samples --out_dir $out_dir --reductionratio $reductionratio --expansionratio $expansionratio

		# We delete the trees files (they take a lot of memory space)
		rm -f ${out_dir}/ts_${repet}_${generationmitoses}.trees

		#python3 sfs_histogram.py --out_dir $out_dir --generationmitoses $generationmitoses --replicatenumber $repet


		# -------------------------------------------------------------------------------------------------------------------------------

		# ------------------ We modify the VCF file created so that it becomes a compatible VCF with various softwares ------------------

		# ---------------------------------------------- and compress it ----------------------------------------------------------------
		
		path="$out_dir/vcf/ts_${repet}_${generationmitoses}.vcf"  # path to the original VCF (created by our call to tree_to_vcf.py)
		if [ $keep_vcf == "false" ]; then  # If there is no need to keep the original VCF (no need to create the .modified files)
				#rm -f "$out_dir/vcf/ts_${repet}_${generationmitoses}.modified.vcf" # we supress the previous .modified file if it existed (to save space)
				#rm -f "$out_dir/vcf/ts_${repet}_${generationmitoses}.modified.vcf.gz" # and the corresponding compressed file
				output=$path # path to the new VCF (slight handmade modifications of the original VCF)
			else # If we need to keep the original VCF file, we save the modifications on a .modified file
				output="$out_dir/vcf/ts_${repet}_${generationmitoses}.modified.vcf" # path to the new VCF (slight handmade modifications of the original VCF)
			fi
		
		# Handmade modifications :
		#awk 'BEGIN {OFS="\t"} { if ($8 == ".") $8 = "DP=1000"; print ; {OFS="\t"}}' $path > t1.vcf  # Transformation of the 8eme (INFO) column from . to DP=1000
		#awk 'BEGIN {OFS="\t"} { if (int($2) >= 1000000)  $1 = 2 ; print ; {OFS="\t"} }' t1.vcf > t2.vcf  # Transformation of the contig name from 1 to 2 in the 1st column (CHROM) if the element in the 2nd column (POS) is greater or equal to 1million
		
		# If we encounter the motif 'contig' we print contig 1 and contig 2, if we don't we print the ligne as usual 
		#awk '/contig/ {print "##contig=<ID=1,length=1000000>" ; print "##contig=<ID=2,length=1000000>"} !/contig/ {print}' t2.vcf > t3.vcf
		# We add the ##INFO ligne between the old ligne number 2 and 3
		#awk 'NR==2{print; print "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"" "Hand modified file" "\">"} NR!=2' t3.vcf > $output

		#rm -f t1.vcf t2.vcf t3.vcf

		# We compress the modified VCF
		#bgzip -f  $output -k #-k is to keep the original .vcf file

		# -------------------------------------------------------------------------------------------------------------------------------

		
		
		done
done


