#!/bin/bash

atoms=$1
ECUT=$2


data_dir="Results_${atoms}"

# 2. Create the folder if it doesn't exist (-p prevents error if it already exists)
if [ ! -d "$data_dir" ]; then
    mkdir -p "$data_dir"
fi

data="./${data_dir}/energy_summary_convergence_kmesh_${atoms}_${ECUT}.txt"

# Check if energy_summary.dat exists and delete it if so
#if [ -f "$data" ]; then
#   rm "$data"
#fi

export ECUT
export atoms 

for ks in {2..12..2} 
do

export ks
file="./${data_dir}/${atoms}$ks"

envsubst < ./template/template.${atoms}.in > ${file}.$ECUT.scf.in


# Run the SCF calculation
pw.x < $file.$ECUT.scf.in > $file.${ECUT}.scf.out

# Extract the total energy and append to the summary file
# (Note: I added quotes around variables to be safe)
ETOT=$(grep '!' $file.${ECUT}.scf.out | cut -d'=' -f2 | cut -d'R' -f1)
echo "$ks $ETOT" >> "$data"

done
