#!/bin/bash
                                                                                     
## Lab : LISN ~ UMR9015 ~ BIOINFO team

#SBATCH --job-name=sfs_representation
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

# Example of working calls :
# sbatch sln_model.sh TYPE="reduction" SWEEP="false"
# sbatch sln_model.sh TYPE="reduction" SWEEP="true"


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


if [ -v REP ]; then
	echo "REP = $REP"
	replicatenumber=$REP  
else 
	replicatenumber=1 
fi

python3 sfs_representation.py --replicatenumber $replicatenumber --type $TYPE