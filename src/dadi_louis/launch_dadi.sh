#!/bin/bash
                                                                                     
## Lab : LISN ~ UMR9015 ~ BIOINFO team

#SBATCH --job-name=dadi_louis
#SBATCH --array=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --account=inferscerevol

# Necessary libraries 
module load bcftools/1.16 
module load dadi/2.2.0 
module load python/3.12 #python/3.10 # test to try and debug smc++
module load mysql-client/8.0.23 
module load psmc/0.6.5 
module load smcpp/1.15.4 
module load vcftools/0.1.16
# no module load available for pyyaml
pip install pyyaml

# ---------------- recup args --------------------------
for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   export "$KEY"="$VALUE"
done
echo "CONFIG_FILE = $CONFIG_FILE"

python3 script_dadi_louis.py --config_file $CONFIG_FILE