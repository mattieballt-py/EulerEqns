import numpy
import matplotlib.pyplot as pl

x = 3
y = 3

h = 0.001
z = x + 2

import numpy as np
import matplotlib.pyplot as plt  

"""
Constants
"""
U = 10 # ms-1
v = 1e-3 # Pas, dynamic viscocity of water, https://www.omnicalculator.com/physics/water-viscosity





"""
Set up of BC's
"""

# want to set y = 0, all t u(0,t)= Ucoswt


"""
Code
"""
# u is a function of y,t so an array with t rows and 
# need a grid in x and then the y is time

# parameters
a = 0  # m  (the x and y boundaries of the bar length along x)
b = 0.5  # m
dx = 0.008  # m  (h, the spatial step)

N = int((b - a) / dx)  # Number of spatial points
x = np.linspace(a, b, N + 1)  # x has length N+1, going from a to b

tend = 3600  # s (simulation time: 1 hour)
dt = 2 # time step
P = int(tend / dt)  # P is Nt, the number of time steps
t = np.linspace(0, tend, P + 1)  # t has length P+1, going from 0 to tend

Ta = 50  # °C temp at x = a for all of t
Tb = 70  # °C temp at x = b for all of t

alpha = 1.172 * 10**-5  # thermal diffusivity constant

# Initialize the domain grid in x-t
Domain = np.zeros((N + 1, P + 1))  # Create an empty domain grid in x and t
Domain[:, 0] = 10  # Initial condition: the entire bar is at 10°C initially
Domain[0, :] = Ta  # Set boundary condition: temp at x = a remains 50°C
Domain[-1, :] = Tb  # Set boundary condition: temp at x = b remains 70°C

# Debugging print statements
print("Shape of Domain:", Domain.shape)

# Stability condition check: required for stability of the numerical solution
stability_factor = alpha * dt / dx**2
print(f"Stability factor = {stability_factor}")
if stability_factor > 0.5:
    raise ValueError("Stability condition violated! Reduce dt or increase dx.")

# Function to compute the temperature evolution
def Temps(Domain, alpha, dt, dx, N, P):
    T = Domain.copy()  # Create a copy of the Domain grid to store updates

    for p in range(1, P + 1):  # Time-stepping loop
        for i in range(1, N):  # Spatial loop (ignores boundary points)
            # Apply the finite difference scheme:
            # T[i, p] = new temperature at spatial point i and time step p
            T[i, p] = ((alpha * dt) / (dx**2)) * (T[i + 1, p - 1] + T[i - 1, p - 1]) + \
                      (1 - 2 * (alpha * dt) / (dx**2)) * T[i, p - 1]

    return T