#!/bin/bash


read -p "What is your atoms name? " atoms
data_dir="Results_${atoms}" 

read -p "Do you want to see the convergence of Ecut? " c_ecut

data_ecut=None
data_kmesh=None
ks=None
Ecut=None

if [ $c_ecut == 1 ]; then
	echo c_ecut
	read -p "What is your k-mesh? " ks
	data_ecut="./${data_dir}/energy_summary_${atoms}${ks}.txt"
fi

read -p "Do you want to see the convergence of kmesh ? " c_kmesh

if [ $c_kmesh == 1 ]; then
	read -p "What is your Ecut? " Ecut
	data_kmesh="./${data_dir}/energy_summary_convergence_kmesh_${atoms}_${Ecut}.txt"
fi

export ks
export data_kmesh
export Ecut
export data_ecut
export data_dir
export data
export c_ecut
export c_kmesh
export atoms

envsubst < ./template/template.plot_deltaE.py > plot_deltaE.py
envsubst < ./template/template.plot_ecut.py > plot_ecut.py

python3 plot_ecut.py
python3 plot_deltaE.py
