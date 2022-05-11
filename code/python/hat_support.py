import numpy as np
from scipy.stats import norm
from numba import jit

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

def expanding_grid(a, b, θ, n):
    grid = np.linspace(0, 1, n)

    grid_e = np.zeros(n)
    for x_i, x in enumerate(grid): 
        grid_e[x_i] = a + (b-a)*x**θ

    return grid_e


