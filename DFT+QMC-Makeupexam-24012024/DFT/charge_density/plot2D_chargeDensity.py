import numpy as np
import matplotlib.pyplot as plt
import matplotlib


atoms = "MgB2"

# Load the data (x, y, f(x, y)) from the file
data = np.loadtxt(f'{atoms}.gnuplot')

# Extract x, y, and f(x, y) values
x = data[:, 0]
y = data[:, 1]
z = data[:, 2]

# Create a grid for the x and y values
nx = int(np.sqrt(len(x)))  # Number of points along one axis (assuming square grid)
ny = nx
x_grid = x.reshape((nx, ny))
y_grid = y.reshape((nx, ny))
z_grid = z.reshape((nx, ny))

# Define saturation limits for z
z_min = 0.0  # Lower bound (e.g., minimum charge density to display)
z_max = 0.5  # Upper bound (saturate above this value)

# Create the saturated colormap contour plot
plt.figure(figsize=(8, 6))
contour = plt.contourf(
    x_grid, y_grid, z_grid, levels=np.linspace(z_min, z_max, 300), cmap='plasma'
)



#contour = plt.contourf(x_grid, y_grid, z_grid, levels=50, cmap='plasma', vmin=z_min, vmax=z_max)
plt.colorbar(contour, label='Charge Density (e/Bohr³)')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'{atoms}')

# Update tick labels to range from -1 to 1
plt.xticks(np.linspace(x_grid.min(), x_grid.max(), 5), np.linspace(-1, 1, 5))
plt.yticks(np.linspace(y_grid.min(), y_grid.max(), 5), np.linspace(-1, 1, 5))
plt.savefig(f'{atoms}_charge_density_vacancy.png', dpi=300, bbox_inches='tight')
plt.show()


