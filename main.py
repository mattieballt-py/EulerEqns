# In this section I am importing all the libraries I will need
import numpy as np
import matplotlib.pyplot as plt  

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

y = np.linspace(0, L, Ny) # for plotting u against the jth space step
t = np.linspace(0, T, Nt) # for plotting u against the nth time step


for n in range(Nt): # going from 0 to Nt-1
    # In this section I am setting the boundary conditions/initial values
    u[0] = np.cos(n * dt)  # Oscillating plate at y = 0
    u[-1] = 0              # Far-field boundary condition at y=L, u=0

    # In this section I am implementing the numerical method
    for j in range(1, Ny - 1):
        u_new[j] = r*u[j + 1] + r*u[j - 1]
    
    u[:] = u_new[:] # Update the solution for the next time step

    if n % 100 == 0:  # Plot every 100 time steps by check if n is divisible by 100
        current_time = round((n * dt),1) # Get the current time and round to 1 decimal point
        label_text = "t = " + str(current_time) # Create a label that shows the time at this step
        plt.plot(y, u, label=label_text) # Plot the current solution with the corresponding label


# Plot results
plt.xlabel('y')
plt.ylabel('u(y, t)')
plt.title("Stokes' Second Problem: Velocity Profile Over Time")
plt.legend()
plt.show()




# In this section I am celebrating
print('CW done: I deserve a good mark')
