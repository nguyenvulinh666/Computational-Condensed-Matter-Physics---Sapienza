
atoms=$1
ks=$2
data_dir="Results_${atoms}${ks}"

export atoms
export data_dir
export ks

envsubst < ./template/template.rho.pp.in > ./charge_density/rho_${atoms}${ks}.pp.in 

pp.x < ./charge_density/rho_$atoms${ks}.pp.in > ./charge_density/rho_$atoms${ks}.pp.out


