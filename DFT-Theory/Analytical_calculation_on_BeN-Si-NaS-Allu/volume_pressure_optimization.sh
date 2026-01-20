#!/bin/bash

atoms=$1

grep vol $atoms.vc_relax.out | cut -d'=' -f2 | awk '{print $1} ' | head -n -1  > $atoms.volume.vc_relax.txt
grep "P=" $atoms.vc_relax.out | cut -d"=" -f2 > $atoms.pressure.vc_relax.txt

paste $atoms.volume.vc_relax.txt $atoms.pressure.vc_relax.txt | column -t >> $atoms.vol_pres.vc_relax.txt
