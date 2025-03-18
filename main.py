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
t = np.arange(0,5,dt)  # x has length N+1, going from a to b
y = np.arange(0,10,dy) # y array with 
N = len(y)
P = len(t)

# Initialize the domain grid in x-t
Domain = np.zeros((N,P))  # Create an empty domain grid in x and t


"""
Set up of BC's
"""
# want to set y = 0, all t u(0,t)= Ucoswt
Domain[:,0] = U * np.cos(w * t)  # °C temp at x = a for all of t
Domain[:,-1] = 0  # ms-1 at last y for all of t = 0

"""
Gauss
"""
# Gauss elimination:
def MyGauss(A,b):
    # number of equations
    n = len(b)
    # eliminate the unknowns, from first to (n-1)th unknown, to form an upper triangular matrix
    for i in range(0,n-1):
        # eliminate the i-th unknown from the (i+1)th row downwards
        for j in range(i+1,n):
            # A(i,i) is the pivot coefficient
            p = A[j,i] / A[i,i]
            # A[j,k] = A[j,k] - p * A[i,k]
            A[j,:] -= p * A[i,:]
            # compute the new element of row j in vector b
            b[j] -= p * b[i]
    # evaluate, by back substitution the solution
    x = np.zeros(n)
    for i in range(n-1,-1,-1):
        # contribution from b (right hand side of the equation)
        x[i] = b[i] / A[i,i]
        for k in range(i+1,n):
            x[i] = x[i] - A[i,k] * x[k] / A[i,i]
    return x

"""
Solving PDE using Crank-Nicholson
"""
def u_cranksolve(Domain, v, dy, dt, N, P):
    r = v * dt / (2 * dy**2)
    A = np.zeros((N-2, N-2))
    b = np.zeros(N-2)
    
    # Construct coefficient matrix A
    for i in range(N-2):
        A[i, i] = 4 # Main diagonal
        if i > 0:
            A[i, i-1] = -r  # Lower diagonal
        if i < N-3:
            A[i, i+1] = -r  # Upper diagonal
    
    # Time-stepping loop
    for k in range(0, P-1):
        # Construct RHS vector b
        for i in range(1, N-1):
            b[i-1] = (1 - 2*r) * Domain[i, k] + r * (Domain[i+1, k] + Domain[i-1, k])
        
        # Apply boundary conditions
        b[0] += r * Domain[0, k]  # Lower BC
        b[-1] += r * Domain[-1, k]  # Upper BC
        
        # Solve using Gauss elimination
        Domain[1:-1, k+1] = MyGauss(A.copy(), b.copy())
    
    return Domain


# solve PDE using function
u = u_cranksolve(Domain, v, dy, dt, N, P)

"""
Scatter Plot of Velocity Profiles at Specific Times
"""
time_indices = [0, int(1/dt), int(5/dt)]  # Indices for t=0s, t=1s, t=5s
time_labels = ["t = 0s", "t = 1s", "t = 5s"]

plt.figure(figsize=(8, 6))
for i, idx in enumerate(time_indices):
    plt.scatter(u[:, idx], y, label=time_labels[i])

plt.xlabel("Velocity (m/s)")
plt.ylabel("Depth (m)")
plt.title("Velocity Profiles at Different Times")
plt.legend()
plt.grid()
plt.show()
