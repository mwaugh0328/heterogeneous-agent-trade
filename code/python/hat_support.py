import numpy as np
import scipy as sp
from scipy.stats import norm
from numba import jit

def expanding_grid(a, b, θ, n):
    grid = np.linspace(0, 1, n)

    grid_e = np.zeros(n)
    for x_i, x in enumerate(grid): 
        grid_e[x_i] = a + (b-a)*x**θ

    return grid_e

def tauchen(n, ρ, σ, m=3, tol=1e-8):
        """
        Inputs: 
        -------------
        NAME        TYPE            DESCR
        n           float64         number of states   
        ρ           float64         AR(1) coefficient    
        σ           float64         Standard deviation of AR(1) innovation
        m           float64         The number of standard deviations to approximate out to
        tol         float64         tolerance for stationary distribution 

        Outputs: 
        -------------
        NAME        TYPE            DESCR
        z           np.array        States of the MArkov chain
        P           np.array        (n x n) transition matrix for Markov chain
        zstar       np.array        stationary distribution 
        """
 
        P = np.empty([n, n])
        z = np.empty(n)
        z_hat = np.empty(n-1)
    
        z_max = m * np.sqrt(σ ** 2 / (1 - (ρ ** 2))) 
        z_min = (-1) * z_max

        # construct n-1 invervals between z_min and z_max
        for i in range(n):
            z[i] = z_min + ((z_max-z_min)/(n-1))*i        

        # construct midpoints
        for i in range(n-1):
            z_hat[i] = (1/2) * (z[i+1]+z[i]) 

        for i in range(n):
            for j in range(n):

                if j == 0:
                    P[i,j] = norm.cdf((z_hat[j]-ρ*z[i])/σ)
                elif j == n-1:
                    P[i,j] = 1 - norm.cdf((z_hat[j-1]-ρ*z[i])/σ)
                else:
                    P[i,j] = norm.cdf((z_hat[j]-ρ*z[i])/σ) - norm.cdf((z_hat[j-1]-ρ*z[i])/σ)

        # stationary distribution: z_star = P*z_star 
        err = 1
        Pm=P
        while err > tol: 
            Pp=Pm@P
            err = np.max(np.abs(Pp-Pm)) 
            Pm=Pp
        zstar = Pp[0,:]

        return z, P, zstar 

def rouwenhorst(n, ρ, σ, μ = 0.0):
    σ_y = σ / (1 - ρ**2)**0.5
    p = (1 + ρ)*0.5
    Θ = np.array([[p, 1 - p], [1 - p, p]])
    ψ = ((n - 1)**0.5) * σ_y
    m = μ / (1 - ρ)
    z, P = _rouwenhorst(p, p, m, ψ, n)

    return z, P 

def _rouwenhorst(p, q, m, Δ, n):
    if n == 2:
        return np.array([m - Δ, m + Δ]),  np.array([[p, 1 - p], [1 - q, q]])
    else:
        _, θ_nm1 = _rouwenhorst(p, q, m, Δ, n - 1)

        θn = (p * np.append(np.append(θ_nm1, np.zeros((n-1, 1)), 1), np.zeros((1, n)), 0) +
             (1 - p) * np.append(np.append(np.zeros((n-1, 1)), θ_nm1, 1), np.zeros((1, n)), 0) +
             q * np.append(np.zeros((1, n)), np.append(np.zeros((n-1, 1)), θ_nm1, 1), 0) +
             (1 - q) * np.append(np.zeros((1, n)), np.append(θ_nm1, np.zeros((n-1, 1)), 1), 0)) 

        θn[1:-1, :] /= 2

        return np.linspace(m - Δ, m + Δ, n), θn

def markov_chain(n_p, n_iid, ρ_p, σ_p, σ_iid, w0=1, tol=1e-8):
    # adapted from Manuel Amador's julia code
    # numbers from KMP
    # ar1 = 0.9695
    # sigmaP = sqrt(0.0384)/(1.2)
    # sigmaIID = sqrt(0.0522)/(1.2)
    # P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)

    # persistent component
    y_vals, P_p = rouwenhorst(n_p, ρ_p, σ_p)
    y_vals = np.exp(y_vals)

    # iid component
    if n_iid > 1:
        x_vals, P_iid = rouwenhorst(n_iid, 0.0, σ_iid)
        x_vals = np.exp(x_vals)
    else:
        P_idd = 1.0
        x_vals = 1.0

    # Initialize arrays
    n_z = n_p * n_iid
    z_vals = np.empty(n_z) 
    P = np.empty((n_z, n_z))

    idx = 0
    for i in range(n_p):
        for j in range(n_iid):
            z_vals[idx] = y_vals[i] * x_vals[j]
            idx += 1

    idx = 0
    for i in range(n_p):
        for j in range(n_iid):
            idx2 = 0
            for k in range(n_p):
                for l in range(n_iid):
                    P[idx,idx2] = P_p[i,k] * P_iid[j,l]
                    idx2 += 1
            idx += 1

    # Find stationary distibution and demean
    # Find eigenvalue closest to 1
    w, v = sp.sparse.linalg.eigs(P.T, k=1, sigma=1)

    # Find stationary distribution
    zstar = np.real(v.flatten())
    zstar /= zstar.sum()    

    # Find ergodic mean
    z_mean = z_vals @ zstar 
    z_vals = (w0/z_mean) * z_vals

    return np.log(z_vals), P