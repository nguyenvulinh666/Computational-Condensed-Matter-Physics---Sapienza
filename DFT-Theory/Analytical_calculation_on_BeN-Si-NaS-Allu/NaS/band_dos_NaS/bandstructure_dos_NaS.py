import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
import matplotlib

atoms="NaS"
Ef =   1.462 # Remember to update the value!
ne = 15
emin = -4
emax = 10
dos_max = 20
num_bands = 10


x_positions = [0.0000, 1.0000, 1.5000, 2.2071, 3.0731]
x_labels = [r'$\Gamma$', 'X', 'W', 'L', r'$\Gamma$']

fig = plt.figure(figsize=(12, 6))
fontsize = 12

matplotlib.rc('xtick', labelsize=13) 
matplotlib.rc('ytick', labelsize=13) 
# Use GridSpec to control width ratios (3 parts left, 1 part right)
gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1]) 
# Set wspace to 0.05 to mimic the tight grouping in your image
plt.subplots_adjust(wspace=0.05) 

# --- LEFT BLOCK: Band Structure ---
ax1 = plt.subplot(gs[0])

bands = np.genfromtxt( f'{atoms}.bands.dat.gnu' )
## Shift all energies by Fermi energy (can be extracted from QE output files)
bands[:,1] = bands[:,1] - Ef
fig = plt.figure(1, [10, 6], dpi=100)

ax1.set_title(f'{atoms} Band structure', fontsize=fontsize+2)
ax1.set_ylabel(r'Energy (eV)', fontsize=fontsize+2)
ax1.set_xlim(0.0,  3.0731)
ax1.set_ylim(emin, emax)
## ## Special k points
ax1.set_xticks(x_positions, x_labels)
## Coordinates can be extracted from the output of bands.x
#ax1.set_yticks(np.linspace(emin, emax, num=8, dtype=float, endpoint=True))
#plt.yticks(np.linspace(-1,0, num=8, dtype=float, endpoint=True))
ax1.grid(axis='x', color='k', linestyle='-')
ax1.locator_params(axis="y", nbins=6)
#Number of bands and total number of k-points + 1 along the chosen path
num_pts = 161

for i in range(num_bands):
    ax1.plot(bands[i*num_pts:i*num_pts+num_pts,0], bands[i*num_pts:i*num_pts+num_pts,1], linestyle = '-', linewidth=2.0)
## Plot Fermi level
ax1.axhline(y = 0.0, color = 'b', linestyle = '--')



# PLot 2
dos = np.genfromtxt(f'{atoms}.dos') 
dos[:,0] = dos[:,0] - Ef
	
ax2 = plt.subplot(gs[1], sharey=ax1) 
	
# Plot DOS (Note: x and y are swapped here to plot vertically)
ax2.plot(dos[:,1]*10, dos[:,0], color = 'k', alpha=0.8, linewidth=1, label=r"D($\epsilon$)")
ax2.plot(dos[:,2], dos[:,0], color = 'r', alpha=0.8, linewidth=1, label=r"N($\epsilon$)")
ax2.axhline(0, linestyle='--', color='b', label=r"EF") 
ax2.vlines(ne, emin, emax, linestyle='--', color='b', label=r"ve")
	
# Add a "Free electron" like curve in red (just for looks)
# ax2.plot(np.sqrt(np.abs(energy_range + 12))/2, energy_range, 'r-', linewidth=1)
# Fermi Level (Blue dashed line)
	
# Styling Right Block
ax2.set_xlabel('DOS (st/eV)', fontsize=fontsize+2)
ax2.set_xlim(0, dos_max)
ax2.set_ylim(emin, emax) 

# Turn off Y-ticks on the left side of this plot (since they are on ax1)
# And move them to the right side
ax2.yaxis.tick_right()
ax2.yaxis.set_label_position("right")
plt.setp(ax2.get_yticklabels(), visible=False)
	
ax2.legend(loc='upper right', bbox_to_anchor=(1.50, 1), borderaxespad=0., ncol=1)
# Final Layout Adjustment
	
plt.savefig(f'{atoms}_Bandstructure_dos.png', transparent=False, dpi=300, bbox_inches='tight')
plt.show()
