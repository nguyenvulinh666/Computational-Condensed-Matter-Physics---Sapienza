### Insstall libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

atoms = "BeN"
num_bands = 10

bands = np.genfromtxt(f'{atoms}.bands.dat.gnu' )
## Shift all energies by Fermi energy (can be extracted from QE output files)
Ef =  9.089 # Remember to update the value!
bands[:,1] = bands[:,1] - Ef
fig = plt.figure(1, [10, 6], dpi=100)

plt.title(f'{atoms} Band structure')
plt.xlim(0.0,  6.0656)
#plt.ylim(-20,20)
#plt.ylim(-1,0.125)
## ## Special k points
plt.xticks([0.0000, 1.0000, 1.5000, 1.8536, 2.9142, 3.7802, 4.3926, 4.7462, 5.4533, 6.0656, 6.0656, 6.4192],[r'$\Gamma$', 'X', 'W', 'K', r'$\Gamma$', 'L', 'U', 'W', 'L', '', 'K|U',  'X']) 
## Coordinates can be extracted from the output of bands.x
#plt.yticks(np.linspace(-20, 20.0, num=8, dtype=float, endpoint=True))
#plt.yticks(np.linspace(-1, 0.125, num=8, dtype=float, endpoint=True))
plt.grid(axis='x', color='k', linestyle='-')
plt.ylabel(r'Energy ($eV$)')

matplotlib.rc('xtick', labelsize=13) 
matplotlib.rc('ytick', labelsize=13) 

#Number of bands and total number of k-points + 1 along the chosen path
num_pts = 161

for i in range(num_bands):
    plt.plot(bands[i*num_pts:i*num_pts+num_pts,0], bands[i*num_pts:i*num_pts+num_pts,1], linestyle = '-', linewidth=2.0)
## Plot Fermi level
plt.axhline(y = 0.0, color = 'b', linestyle = '--')
plt.savefig(f'{atoms}.png', transparent=False, dpi=300, bbox_inches='tight')
plt.show()
