#!/bin/bash

atoms=$1
ks=$2
nat=$3
ECUT=$4

prefix=${atoms}${ks}

data_dir="bands_dos_${prefix}"

# 2. Create the folder if it doesn't exist (-p prevents error if it already exists)
if [ ! -d "$data_dir" ]; then
    mkdir -p "$data_dir"
fi

data="./${data_dir}/energy_summary_${atoms}${ks}.txt"


export ks
export nat
export atoms
export ECUT 
export prefix

pw.x < ./Results_${prefix}/${prefix}.${ECUT}.scf.in > ./Results_${prefix}/${prefix}.${ECUT}.scf.out

ks=$(($ks+8))

envsubst < ./template/template.${atoms}.nscf.in > ./${data_dir}/${prefix}.nscf.in
pw.x < ./${data_dir}/${prefix}.nscf.in > ./${data_dir}/${prefix}.nscf.out


#ECUT=$(($ECUT * 2))
envsubst < ./template/template.${atoms}.bands.in > ./${data_dir}/${prefix}.bands.in
#envsubst < ./template/template.bands_${atoms}.in > ./${data_dir}/bands_${prefix}.in

pw.x < ./${data_dir}/${prefix}.bands.in > ./${data_dir}/${prefix}.bands.out

# Create bands file
cat > ./${data_dir}/bands_${prefix}.in << EOF
 &BANDS                                                                
 prefix = '${prefix}',                                                                  
 outdir = './tmp/',                                                                     
 filband = '${prefix}.bands.dat',                                                
 lsym = .true.                                                                          
 &end 
EOF

bands.x < ./${data_dir}/bands_${prefix}.in > ./${data_dir}/bands_${prefix}.out

# Create dos file
cat > ./${data_dir}/dos_$prefix.in << EOF
&dos
prefix = '${prefix}',
outdir = './tmp',
DeltaE = 0.01
&end
EOF

dos.x < ./${data_dir}/dos_$prefix.in > ./${data_dir}/dos_$prefix.out

mv $prefix.dos ./${data_dir}/$prefix.dos
mv $prefix.bands.dat.rap ./${data_dir}/$prefix.bands.dat.rap
mv $prefix.bands.dat.gnu ./${data_dir}/$prefix.bands.dat.gnu
mv $prefix.bands.dat ./${data_dir}/$prefix.bands.dat


