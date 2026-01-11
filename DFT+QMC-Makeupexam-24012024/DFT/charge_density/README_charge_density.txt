# Command, you should change all the name of lattice inside the file, so it can run probably
pp.x < rho_MgB2.pp.in > rho.pp.out # Plot charge density
pp.x < charge_density_MgB2.in > charge_density_MgB2.out # Plot charge density on the cut
python3 plot2D_chargeDensity.py 


