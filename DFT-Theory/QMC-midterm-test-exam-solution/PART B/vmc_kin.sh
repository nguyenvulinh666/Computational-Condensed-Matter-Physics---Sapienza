mkdir "$1"
cd "$1"
tar -xf ../exam_bin.tar.gz
cp ../vmc_kin.sh ./
#gfortran setup.f -O3 -w -o setup.x
#gfortran qmc.f -O3 -w -std=legacy -o qmc.x
#gfortran statfor.f -O3 -w -o statfor.x
#gfortran statforw.f -O3 -w -o statforw.x
rs=20     # r_s paremeter of the wavefunction
nup=21    # number of spin up electrons
ndw=0    # number of spin down electrons
twb=1     # whether to use 2-body correlation functions
thb=0     # whether to use 3-body correlation functions
bkf=0     # whether to use backflow correlation functions
blocks=200
steps=1000
timestep=4
cat << EOF > input
${1}
${rs}
${nup}
${ndw}
${twb}
${thb}
${bkf}
EOF
./setup.x < input > "$1".setup
if [ "$thb" -eq "1" -o "$bkf" -eq "1" ]
then
echo "OPTIMIZING"
# 500 is the number of steps used for every optimization, and 10 is the frequency of saving configuration
cat << EOF > "$1".in
vmc 100 50 10 10
optimize 500 ${1}.b ${1}.u3
EOF
info=-1
while [ "$info" -eq "-1" ]
do
./qmc.x < runid | tee -a "$1".optmz.out >/dev/null
info=$(cat "$1".optmz.out | grep info | tail -1 | awk '{ print $5 }')
echo $info
done
mv "$1".in "$1".optmz.in
fi
echo "SIMULATING"
cat << EOF > "$1".in
vmc ${blocks} ${steps} ${timestep}
EOF
./qmc.x < runid > "$1".out
grep elocal "$1".out > "$1".elocal
grep ekin "$1".out > "$1".ekin
grep "acc.rate" "$1".out > "$1".accrate
cat "$1".elocal | ./statfor.x > "$1".elocal.stat
cat "$1".ekin | ./statfor.x > "$1".ekin.stat
cat "$1".accrate | ./statfor.x > "$1".accrate.stat
echo "ELOCAL"
cat "$1".elocal.stat
echo "EKIN"
cat "$1".ekin.stat
echo "ACC.RATE"
cat "$1".accrate.stat
gnuplot <<-EOFMarker
    set term png
    set output '${1}.elocal.png'
    plot '${1}.elocal' using 0:1 with linespoints
EOFMarker
cd ..
tar -czf "$1".tar.gz ./"$1"/*
rm -rf "$1"
