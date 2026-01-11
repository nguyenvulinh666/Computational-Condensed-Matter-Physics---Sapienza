#!/bin/bash

atoms=$1
ks=$2


data_dir="Results_${atoms}"

# 2. Create the folder if it doesn't exist (-p prevents error if it already exists)
if [ ! -d "$data_dir" ]; then
    mkdir -p "$data_dir"
fi

file="./${data_dir}/${atoms}$ks"
data="./${data_dir}/energy_summary_${atoms}${ks}.txt"

# Check if energy_summary.dat exists and delete it if so
if [ -f "$data" ]; then
   rm "$data"
fi

export ks
export atoms 

for ECUT in 20 40 60 80 100 120 
do

export ECUT

envsubst < ./template/template.${atoms}.in > ${file}.$ECUT.scf.in


# Run the SCF calculation
pw.x < $file.$ECUT.scf.in > $file.${ECUT}.scf.out

# Extract the total energy and append to the summary file
# (Note: I added quotes around variables to be safe)
ETOT=$(grep '!' $file.${ECUT}.scf.out | cut -d'=' -f2 | cut -d'R' -f1)
echo "$ECUT $ETOT" >> "$data"

done
