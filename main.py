import numpy as np
import matplotlib.pyplot as plt  

"""
Constants
"""
U = 10 # ms-1, amplitude of plate oscillations
v = 1e-3 # Pas, dynamic viscocity of water, https://www.omnicalculator.com/physics/water-viscosity
w = 10 # Hz, frequency of oscillations
tend = 3600  # s (simulation time: 1 hour)

# parameters
dt = 0.05 # s, step size in t grid
dy = np.sqrt(dt*v) # setting r = 1, so dy changes to allow this
print("dy now",dy)
# Arrays
t = np.arange(0,5,dt)  # x has length N+1, going from a to b
y = np.arange(0,10,dy) # y array with 
N = len(t)
P = len(y)
# Initialize the domain grid in x-t
Domain = np.zeros((len(y), len(t)))  # Create an empty domain grid in x and t


"""
Set up of BC's
"""
# want to set y = 0, all t u(0,t)= Ucoswt
Domain[0,:] = U * np.cos(w * t)  # °C temp at x = a for all of t
Domain[-1,:] = 0  # ms-1 at last y for all of t = 0

"""
Stability Condition Check
"""
r = 1
# Stability condition check: required for stability of the numerical solution
stability_factor = r
print("Stability factor =", stability_factor)
if stability_factor > 0.5:
    raise ValueError("Stability condition violated! Reduce dt or increase dx.")

"""
Solving PDE 
"""
# u is a function of y,t so an array with t rows and 
def u(Domain, v, dy, dt, N, P):
    for k in range(0, len(t) + 1):
        print("hi")
    return


print("t",t,"y",y)

Domain = u
# Debugging print statements
print("Shape of Domain:", u)


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
    
T = Temps(Domain, alpha, dt, dy, len(t), len(y))
plt.figure(figsize=(8, 5))


"""
Plot Contour of Velocity Evolution (t on x-axis, y on y-axis)
"""
plt.figure(figsize=(8, 6))

# Create a contour plot (alternative to imshow)
T, Y = np.meshgrid(t, y)  # Create meshgrid for contour plot
plt.contourf(T, Y, U_t.T, levels=20, cmap="coolwarm")

plt.colorbar(label="Velocity (m/s)")
plt.xlabel("Time (s)")
plt.ylabel("Depth (m)")
plt.title("Velocity Evolution in Stokes' Second Problem")
plt.show()




# Gauss elimination:
def MyGauss(A,b):
    
    # number of equations
    n = len(b)
    
    # eliminate the unknowns, from first to (n-1)th unknown, to form an upper triangular matrix
    for i in range(0,n-1):
        # eliminate the i-th unknown from the (i+1)th row downwards
        # i.e. set the zeros in column i.
        for j in range(i+1,n):
            # eliminate on row j

            # A(i,i) is the pivot coefficient
            p = A[j,i] / A[i,i]
        
            # compute the new elements of row j in matrix A
            # use slicing
            #A[j,:] = A[j,:] - p * A[i,:]
            # or, alternatively, loop for every cell of row j
            #for k in range(i,n):
            #    A[j,k] = A[j,k] - p * A[i,k]
            A[j,:] = A[j,:] - p * A[i,:]

            # compute the new element of row j in vector b
            b[j] = b[j] - p * b[i]
    
    
    # evauate, by back substitution the solution
    # start from the last unknown and go upward till the first unknown
    x = np.zeros(n)
    for i in range(n-1,-1,-1):
        # contribution from b (right hand side of the equation)
        x[i] = b[i] / A[i,i]
        # contribution from the other (already evaluated) unknowns
        # (within the left hand side of the equation)
        for k in range(i+1,n):
            x[i] = x[i] - A[i,k] * x[k] / A[i,i]

    return x

