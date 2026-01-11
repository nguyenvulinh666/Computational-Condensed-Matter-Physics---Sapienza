#!/bin/bash
atoms=$1
ks=$2
prefix=${atoms}${ks}

export atoms 
export prefix

envsubst < ./template/template.bandstructure_dos.py > bandstructure_dos_${prefix}.py

python3  bandstructure_dos_${prefix}.py

