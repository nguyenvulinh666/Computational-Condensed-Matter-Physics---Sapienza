mkdir "$1"
cd "$1"
tar -xf ../exam_bin.tar.gz
cp ../dmc_kin.sh ./
#gfortran setup.f -O3 -w -o setup.x
#gfortran qmc.f -O3 -w -std=legacy -o qmc.x
#gfortran statfor.f -O3 -w -o statfor.x
#gfortran statforw.f -O3 -w -o statforw.x
rs=10     # r_s paremeter of the wavefunction
nup=21    # number of spin up electrons
ndw=0    # number of spin down electrons
twb=1     # whether to use 2-body correlation functions
thb=1     # whether to use 3-body correlation functions
bkf=1     # whether to use backflow correlation functions
blocks=100
steps=200
timestep=0.25
energy0=-0.170
walkers=40
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
optimize 500 ${1}.b ${1}.u3z
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
# 400 = walkers * frequency (last number)
cat << EOF > "$1".in
vmc 1 400 10 10
dmc ${blocks} ${steps} ${timestep} ${walkers} ${energy0}
EOF
./qmc.x < runid > "$1".out
grep elocal "$1".out > "$1".elocal
cat "$1".elocal | tail -450 | ./statforw.x > "$1".elocal.stat
grep ekin "$1".out > "$1".ekin
cat "$1".ekin | tail -450 | ./statforw.x > "$1".ekin.stat
echo "ELOCAL"
cat "$1".elocal.stat
echo "ELKIN"
cat "$1".ekin.stat
cd ..
tar -czf "$1".tar.gz ./"$1"/*
rm -rf "$1"
