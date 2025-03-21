# In this section I am importing all the libraries I will need
import numpy as np
import matplotlib.pyplot as plt  
from mpl_toolkits.mplot3d import Axes3D

# hello, we have done the two methods with the template sub headings as below: enjoy!
"""
Implicit Crank-Nicolson Method
"""
# In this section I am setting the domain of solution and the discretised grid

Ny = 200  # Number of spatial points
Nt = 2000  # Number of time steps
L = 10  # Domain size
T = 40  # Total time
h = L / (Ny - 1)  # Spatial step
k = h**2  # Set time step such that r = 1 
r = k / (2 * h**2)  # Crank-Nicolson parameter 

# In this section I am defining arrays I would need 

A = np.zeros((Ny, Ny))
B = np.zeros((Ny, Ny))

for i in range(1, Ny-1):
    A[i, i-1] = -r
    A[i, i] = 1 + 2*r
    A[i, i+1] = -r

    B[i, i-1] = r
    B[i, i] = 1 - 2*r
    B[i, i+1] = r

# In this section I am setting the boundary conditions/initial values

A[0, 0] = 1 # Boundary conditions: Dirichlet at y=0 (oscillating plate)
B[0, 0] = 1

A[-1, -1] = 1 # Neumann boundary at y -> ∞ (no gradient at far boundary)
B[-1, -1] = 1

u = np.zeros(Ny)# Initial condition: fluid at rest
y_values = np.linspace(0, L, Ny)

# In this section I am implementing the numerical method

u_history = np.zeros((Nt, Ny))

for n in range(Nt): # Time-stepping
    u[0] = np.cos(n * k)  # Apply oscillating boundary condition
    b = B @ u  # Compute RHS
    u = np.linalg.solve(A, b)  # Solve system
    u_history[n, :] = u.copy()  # Store solution for animation


# In this section I am showing the results
plt.figure(1)
fig, ax = plt.subplots(figsize=(6, 5))
line, = ax.plot(y_values, u_history[0], 'b-', lw=2)
ax.set_xlim(0, L)
ax.set_ylim(-1, 1)
ax.set_xlabel("non-dimensional height y/δ")
ax.set_ylabel("non-dimensional Velocity u/U")
title = ax.set_title("Velocity Evolution at ωt = 0")

# --- Manual Animation Loop (Works in Spyder and VS Code) ---
skip = 3 # Set the 'frame rate'

for frame in range(0, Nt, skip):
    line.set_ydata(u_history[frame])  # Update velocity profile
    title.set_text(f"Velocity Profile at ωt = {frame * k:.2f}")  # Update title
    plt.draw()  # Force update
    plt.pause(0.01)  # Allow time for GUI to update


plt.ioff()  # Keep plot open and turn off interactive mode (Spyder)
plt.show()

"""
Explicit Finite Difference Method
"""
# In this section I am setting the domain of solution and the discretised grid

L = 10        # Length of the domain
T = 40        # Total time to simulate
Ny = 200       # Number of spatial points
Nt = 2000    # Number of time steps
dy = L / (Ny - 1)
dt = 0.5*dy**2 # Stability parameter (must satisfy r <= 0.5 for stability)
r = 0.5

# In this section I am defining arrays I would need 

u = np.zeros(Ny) # array of U values at (n) that will be overwritten with updated u values (memory efficient and faster)
u_new = np.zeros(Ny) # array holding new (n + 1) values of U
u_hstry = np.zeros((Nt,Ny)) # 2D array for storing history of U values for plotting surface plot

y = np.linspace(0, L, Ny) # for plotting u against the jth space step
t = np.linspace(0, T, Nt) # for plotting u against the nth time step

# Plot results no.2
plt.figure(2) # to make distinct from the implicit crank nicolson plot
plt.figure(figsize=(6,5))
plt.xlabel('y/δ')
plt.ylabel('u(y, t)/U')
plt.title("Stokes Second Problem solved with Explicit Finite Difference Method")

for n in range(Nt): # going from 0 to Nt-1
    # In this section I am setting the boundary conditions/initial values
    u[0] = np.cos(n * dt)  # Oscillating plate at y = 0
    u[-1] = 0              # Far-field boundary condition at y=L, u=0

    # In this section I am implementing the numerical method
    for j in range(1, Ny - 1):
        u_new[j] = r*u[j + 1] + r*u[j - 1]
    
    u[:] = u_new[:] # Update the solution for the next time step
    u_hstry[n, :] = u[:] # save solution at each time step for plotting
    if n % 100 == 0:  # Plot every 100 time steps by check if n is divisible by 100
        current_time = round((n * dt),1) # Get the current time and round to 1 decimal point
        label_text = "t = " + str(current_time) # Create a label that shows the time at this step
        plt.plot(y, u, label=label_text) # Plot the current solution with the corresponding label



# In this section I am showing the results
plt.legend()
plt.show()

plt.figure(3)
Y, T_grid = np.meshgrid(y, t)  # Create 2D grids for y and t for plotting surface

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

surf = ax.plot_surface(Y, T_grid, u_hstry, cmap='viridis', edgecolor='none') # Plot the surface
ax.set_xlabel('Non-dimensional position (y/δ)')
ax.set_ylabel('Non-Dimensional Time (t*)')
ax.set_zlabel('Non-dimensional velocity u(y, t)/U')
ax.set_title('3D Surface Plot of Diffusion')
cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
cbar.set_label(r'$u/U$')  # the $ needed as it uses Latex formating

plt.show()

# In this section we are celebrating
print('CW done: We deserve a good mark')