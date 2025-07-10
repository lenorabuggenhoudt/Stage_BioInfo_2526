#!/bin/bash
                                                                                     
## Lab : LISN ~ UMR9015 ~ BIOINFO team

#SBATCH --job-name=comparaison
#SBATCH --array=1
#SBATCH --mem-per-cpu=32GB
#SBATCH --account=inferscerevol

# Necessary libraries

module load python/3.12 #python/3.10 # test to try and debug smc++
# -------------------------------------------------------
# ---------------- recup args --------------------------
for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   export "$KEY"="$VALUE"
done
echo "CONFIG_FILE = $CONFIG_FILE"
echo "DADI_LOUIS_FILE = $DADI_FILE"
echo "STAIRWAY_FILE = $STAIRWAY_FILE"
echo "Nu theo = $Nu"
echo "T theo = $T"
echo "SIMU TYPE = $SIMU_TYPE" # modele_selection ex : expansion_neutral


python3 data_extraction_boxplot_test.py --dadi_louis $DADI_FILE --stairway_summary $STAIRWAY_FILE --Nu_theo $Nu --T_theo $T --simu_type $SIMU_TYPE