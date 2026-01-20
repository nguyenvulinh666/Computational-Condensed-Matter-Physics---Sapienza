#!/bin/bash

atoms=$1

grep vol $atoms.vc_relax.2.out | cut -d'=' -f2 | awk '{print $1} ' | head -n -1  > $atoms.volume.vc_relax.2.txt
grep "P=" $atoms.vc_relax.2.out | cut -d"=" -f2 > $atoms.pressure.vc_relax.2.txt

paste $atoms.volume.vc_relax.2.txt $atoms.pressure.vc_relax.2.txt | column -t > $atoms.vol_pres.vc_relax.2.txt
