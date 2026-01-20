import matplotlib.pyplot as plt
import numpy as np
import matplotlib

fontsize = 12
volume = []
pressure = []

atoms = "si.cd"

with open(f"{atoms}.vol_pres.vc_relax.txt", 'r') as file:
        #next(file)  # Skip the header line
        for line in file:
            e, t = line.split()
            volume.append(float(e))
            pressure.append(float(t))


steps = np.linspace(1, len(volume), len(volume))


# 2. Create the Figure and the Primary Axis (Left - Volume)

fig, ax1 = plt.subplots(figsize=(10, 6))




# --- Plot Volume on Left Axis ---

color_vol = 'blue'

ax1.set_xlabel('Relaxation Step', fontsize=fontsize+2)

ax1.set_ylabel('Volume (B$^3$)', color=color_vol, fontsize=fontsize+2) # B^3 likely stands for Bohr^3

lns1 = ax1.plot(steps, volume, color=color_vol, marker='o', linestyle='-', label='Volume (B$^3$)')

ax1.tick_params(axis='y', labelcolor=color_vol)

ax1.grid(True, axis='y') # Add horizontal grid lines

ax1.locator_params(axis="y", nbins=6)
# --- Plot Pressure on Right Axis ---
matplotlib.rc('xtick', labelsize=13) 
matplotlib.rc('ytick', labelsize=13) 

# twinx() creates a second axes that shares the same x-axis
ax2 = ax1.twinx()  


color_press = 'red'

ax2.set_ylabel('Pressure (kbar)', color=color_press, fontsize=14)

lns2 = ax2.plot(steps, pressure, color=color_press, marker='s', linestyle='-', label='Pressure (kbar)')

ax2.tick_params(axis='y', labelcolor=color_press)
ax2.locator_params(axis="y", nbins=6)


# 3. Add Legends

# We add them separately to match the positions in your image

lns = lns1 + lns2
labs = [l.get_label() for l in lns]

ax1.legend(lns, labs, loc='upper center', bbox_to_anchor=(0.5, 1.10), ncol=2, fancybox=True, shadow=False, fontsize=fontsize)


# 4. Title and Layout

# plt.title('Volume and Pressure vs Relative Step', fontsize=16)

plt.tight_layout()

plt.savefig(f'{atoms}_volume_and_pressure_dos.png', transparent=False, dpi=300, bbox_inches='tight')
plt.show()
