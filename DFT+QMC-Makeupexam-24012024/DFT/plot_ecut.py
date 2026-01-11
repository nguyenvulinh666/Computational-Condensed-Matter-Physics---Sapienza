import matplotlib.pyplot as plt
import matplotlib
# Read the data from energy_ecut.dat

c_ecut = 0
c_kmesh = 0

if 1 == 1:
    ecut, etot = [], []
    with open("./Results_MgB2/energy_summary_MgB28.txt", 'r') as file:
        #next(file)  # Skip the header line
        for line in file:
            e, t = line.split()
            ecut.append(float(e))
            etot.append(float(t))

    # Plotting
    plt.figure(figsize=(10, 6))
    matplotlib.rc('xtick', labelsize=13) 
    matplotlib.rc('ytick', labelsize=13) 
    plt.locator_params(axis="both", nbins=6)
    plt.plot(ecut, etot, marker='o', linestyle='-', color='b')
    plt.title('Total Energy (Ryd) vs Ecut (Ryd)', fontsize=18)
    plt.xlabel('Ecut (Ryd)', fontsize=16)
    plt.ylabel('Etot (Ryd)', fontsize=16)
    plt.grid(True)
    plt.xlim(min(ecut) - 10, max(ecut) + 10)  # Adjust x-axis range
    plt.ylim(min(etot) - 0.2, max(etot) + 0.2)  # Adjust y-axis range
    plt.savefig('./Results_MgB2/etot_ecut_MgB28.png')
    plt.show()
    
if 1 == 1: 
    ecut, etot = [], []
    with open("./Results_MgB2/energy_summary_convergence_kmesh_MgB2_90.txt", 'r') as file:
        #next(file)  # Skip the header line
        for line in file:
            e, t = line.split()
            ecut.append(float(e))
            etot.append(float(t))

    # Plotting
    plt.figure(figsize=(10, 6))
    matplotlib.rc('xtick', labelsize=13) 
    matplotlib.rc('ytick', labelsize=13) 
    plt.locator_params(axis="both", nbins=6)
    plt.plot(ecut, etot, marker='o', linestyle='-', color='b')
    plt.title('Total Energy (Ryd) vs Ecut (Ryd)', fontsize=18)
    plt.xlabel('kpt', fontsize=16)
    plt.ylabel('Etot (Ryd)', fontsize=16)
    plt.grid(True)
#    plt.xlim(min(ecut) - 10, max(ecut) + 10)  # Adjust x-axis range
    plt.ylim(min(etot) - 0.2, max(etot) + 0.2)  # Adjust y-axis range
    plt.savefig('./Results_MgB2/etot_kmesh_MgB2_90.png')
    plt.show()
