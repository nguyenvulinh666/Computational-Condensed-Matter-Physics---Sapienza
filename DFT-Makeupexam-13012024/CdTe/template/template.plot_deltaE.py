import matplotlib.pyplot as plt
import matplotlib
import numpy as np

threshold = 20
c_ecut = 0
c_kmesh = 0


# Rydberg → meV conversion factor
RY_TO_MEV = 13605.693
if $c_ecut == 1:	
	# Read the data from energy_ecut.dat
	ecut, etot = [], []
	with open("${data_ecut}", 'r') as file:
	    #next(file)  # skip header
	    for line in file:
	        e, t = line.split()
	        ecut.append(float(e))
	        etot.append(float(t) * RY_TO_MEV)     # convert to meV
	
	# Compute differences: E(i) - E(i+1), in meV
	delta_E = np.array([abs(etot[i] - etot[i+1]) for i in range(len(etot)-1)])
	
	# x-axis at midpoints between consecutive cutoffs
	ecut_mid = [(ecut[i] + ecut[i+1]) / 2 for i in range(len(ecut)-1)]
	
	# Plot
	plt.figure(figsize=(10, 6))
	matplotlib.rc('xtick', labelsize=13) 
	matplotlib.rc('ytick', labelsize=13) 
	plt.locator_params(axis="both", nbins=6)
	
	plt.plot(ecut_mid, delta_E, marker='o', linestyle='-')
	plt.xlabel('Ecut (Ryd)', fontsize=15)
	plt.ylabel(r'|$\Delta E$| (meV)', fontsize=15)
	plt.grid(True)
	
	# Limit y-axis to 500 meV
	plt.ylim(-5, max(delta_E)+(max(delta_E)-min(delta_E))*0.3)
	plt.xlim(min(ecut_mid)-10, max(ecut_mid)+10)
	
	plt.axhline(y = threshold, color ="red", linestyle ="--", label=f"Energy threshold {threshold} meV")
	plt.legend(fontsize=13)
	plt.savefig('./${data_dir}/deltaE_mev_Ecut_${atoms}${ks}.png')
	plt.show()

if $c_kmesh == 1:	
	# Read the data from energy_ecut.dat
	ecut, etot = [], []
	with open("${data_kmesh}", 'r') as file:
	    #next(file)  # skip header
	    for line in file:
	        e, t = line.split()
	        ecut.append(float(e))
	        etot.append(float(t) * RY_TO_MEV)     # convert to meV
	
	# Compute differences: E(i) - E(i+1), in meV
	delta_E = np.array([abs(etot[i] - etot[i+1]) for i in range(len(etot)-1)])
	
	# x-axis at midpoints between consecutive cutoffs
	ecut_mid = [(ecut[i] + ecut[i+1]) / 2 for i in range(len(ecut)-1)]
	
	# Plot
	plt.figure(figsize=(10, 6))
	matplotlib.rc('xtick', labelsize=13) 
	matplotlib.rc('ytick', labelsize=13) 
	plt.locator_params(axis="both", nbins=6)
	
	plt.plot(ecut_mid, delta_E, marker='o', linestyle='-')
	plt.xlabel('kpt', fontsize=15)
	plt.ylabel(r'|$\Delta E$| (meV)', fontsize=15)
	plt.grid(True)
	plt.axhline(y = threshold, color ="red", linestyle ="--", label=f"Energy threshold {threshold} meV")
	plt.legend(fontsize=13)
	# Limit y-axis to 500 meV
	plt.ylim(-5, max(delta_E)+(max(delta_E)-min(delta_E))*0.3)
#	plt.xlim(min(ecut_mid)-10, max(ecut_mid)+10)
	plt.savefig('./${data_dir}/deltaE_mev_kmesh_${atoms}${Ecut}.png')
	plt.show()
