import numpy as np
import scipy as sp
from hat_support import *

class HATmodel:

    def __init__(self, β, σ, α, φ, δ, σ_ε, N, P, L):
        self.β, self.σ, self.α, self.φ = β, σ, α, φ
        self.δ, self.σ_ε, self.N, self.P = δ, σ_ε, N, P
        self.L = L

    def utility(self, c):
        σ = self.σ
        return ((1-σ)**(-1))*(c**(1-σ))

    def u_c(self, c):
        σ = self.σ        
        return c**(-σ)

    def u_c_inv(self, x):
        σ = self.σ        
        
        return x**(-1/σ) 

    def find_w(self, R, TFP):
        α, δ = self.α, self.δ
        return TFP*(1-α)*((R-1+δ)/(TFP*α))**(α/(α-1)) 

    def find_K(self, R, TFP, N):
        α, δ = self.α, self.δ
        return N*((R-1+δ)/(TFP*α))**(1/(α-1)) 

    def find_Y(self, K, TFP, N):
        α = self.α
        return TFP*(K**α)*(N**(1-α)) 

    def find_X(self, K):
        δ = self.δ 
        return δ*K 

    def coleman(self, agrid, zgrid, x, R, w, p):
        u_c, u_c_inv = self.u_c, self.u_c_inv
        β, φ, P, N = self.β, self.φ, self.P, self.N 
        find_value_function = self.find_value_function
        find_probabilities = self.find_probabilities

        # Unpack policy and value function
        gc = x[0:agrid.shape[0], :, :]
        v = x[agrid.shape[0]:, :, :]
 
        # Compute choice probabilities
        π = find_probabilities(v)

        # Integrate over ϵ
        Ev1 = np.zeros((agrid.shape[0], zgrid.shape[0]))
        for a_i, a in enumerate(agrid):
            for z_i, z in enumerate(zgrid):
                Ev1[a_i, z_i] = π[a_i, z_i,:] @ (u_c(gc[a_i, z_i,:])/p)                            

        # Integrate over z
        Ev = Ev1 @ P.T

        gc1 = np.empty_like(gc)
        ga1 = np.empty_like(gc)
        for i in range(N):
           
            # Build endogenous grid assuming interior solution
            endo_c = u_c_inv(p[i]*β*R*Ev)
            endo_a = (1/R)*(p[i]*endo_c+agrid[:,np.newaxis]-w*zgrid)

            # Interpolate to find new policy function
            ga1[:,:,i] = np.array([np.interp(agrid, endo_a[:, z_i], agrid) for z_i in range(len(zgrid))]).T

        # Back out consumption policy
        cih = zgrid*w + R*agrid[:,np.newaxis]
        gc1 = (cih[:,:,np.newaxis] - ga1)/p        

        v1 = find_value_function(agrid, zgrid, v, gc1, ga1)
        x1 = np.concatenate([gc1, v1])

        return x1 

    def find_value_function(self, agrid, zgrid, v, gc, ga):
        find_log_sum, utility = self.find_log_sum, self.utility
        P, N = self.P, self.N
        β = self.β

        for N_i in range(N):
            for z_i, z in enumerate(zgrid):
                for a_i, a in enumerate(agrid): 

                    # Compute Expectations 
                    Ev = 0
                    if ga[a_i, z_i, N_i] <= agrid[0]:
                        Ev += P[z_i, :] @ find_log_sum(v[0, :, :], zgrid)
                    elif ga[a_i, z_i, N_i] >= agrid[-1]:
                        Ev += P[z_i, :] @ find_log_sum(v[agrid.shape[0]-1, :, :], zgrid)
                    else: 
                        idx = np.where(ga[a_i, z_i, N_i]>=agrid)[0][-1]
                        q = ((ga[a_i, z_i, N_i] - agrid[idx]) / (agrid[idx+1] - agrid[idx]))
                        Ev += (1-q) * P[z_i, :] @ find_log_sum(v[idx, :, :], zgrid)
                        Ev += q * P[z_i, :] @ find_log_sum(v[idx+1, :, :], zgrid)

                    v[a_i, z_i, N_i] = utility(gc[a_i, z_i, N_i]) + β*Ev

        return v

    def find_policies(self, agrid, zgrid, R, w, p):
        u_c, u_c_inv = self.u_c, self.u_c_inv
        coleman, find_probabilities = self.coleman, self.find_probabilities      
        φ, β, N = self.φ, self.β, self.N 
        fixed_point_iteration = self.fixed_point_iteration

        # Initial guess for consumption policy
        gc0_vec = np.linspace(0.1, 3, agrid.shape[0])
        gc0 = np.repeat(np.repeat(gc0_vec[:,np.newaxis], zgrid.shape[0], axis=1)[:,:,np.newaxis], N, axis=2) 
        v0 = (-1)*np.ones((agrid.shape[0], zgrid.shape[0], N))/(1-β)
        x0 = np.concatenate([gc0, v0])

        # Find fixed point of coleman operator
        fun = lambda x: coleman(agrid, zgrid, x, R, w, p)     
        x = fixed_point_iteration(fun, x0, agrid)

        gc = x[0:agrid.shape[0], :, :]
        v = x[agrid.shape[0]:, :, :]

        # Find labor and asset policies
        cih = zgrid*w + R*agrid[:,np.newaxis]
        ga = np.maximum(cih[:,:,np.newaxis] - gc*p, φ)       

        # Find probabilities
        π = find_probabilities(v)

        return gc, ga, π

    def fixed_point_iteration(self, fun, x0, agrid, max_iter = 1000, error=1e-8):

        x = x0
        gc = x[0:agrid.shape[0], :, :]
        for i in range(max_iter):
            x1 = fun(x)
            gc1 = x1[0:agrid.shape[0], :, :]

            if np.abs(np.max(gc1-gc)) < error:
                return x1
            else:
                x = x1
                gc = gc1

        print('Fixed point iteration did not converse!')

    def find_probabilities(self, v):
        σ_ε = self.σ_ε

        v_temp = (v.T - np.max(v, axis=2).T).T
        return (np.exp(v_temp/σ_ε).T/np.exp(v_temp/σ_ε).sum(2).T).T

    def find_log_sum(self, vj, zgrid):
        σ_ε = self.σ_ε

        vj_max = np.max(vj, axis=1)
        log_sum = np.array([σ_ε*np.log(np.sum(np.exp((vj[z_i, :]-vj_max[z_i])/σ_ε)))+vj_max[z_i] for z_i in range(len(zgrid))])

        return log_sum

    def find_transition_matrix(self, agrid, zgrid, π, ga):
        """
        Transition matrix for the vector s = [(a1,z1),...,(aM,z1),(a1,z2),...,(aM,z2),...(aM,zN)]'        
        """
        P, N = self.P, self.N

        # Find the transition matrix
        T = np.zeros((agrid.size*zgrid.size, agrid.size*zgrid.size))
        
        for N_i in range(N):
            for z_i, z in enumerate(zgrid):
                for a_i, a in enumerate(agrid):
                    s_i = a_i + z_i*agrid.size

                    if ga[a_i, z_i, N_i] <= agrid[0]:
                        ap_i = 0
                        q = 0
                    elif ga[a_i, z_i, N_i] >= agrid[-1]:
                        ap_i = agrid.shape[0]-2 
                        q = 1
                    else:
                        ap_i = np.where(ga[a_i, z_i, N_i]>=agrid)[0][-1]
                        q = ((ga[a_i, z_i, N_i] - agrid[ap_i]) / (agrid[ap_i+1] - agrid[ap_i]))

                    for zp_i, zp in enumerate(zgrid):
                        sp_i = ap_i + zp_i*agrid.size
                    
                        T[s_i, sp_i] += (1-q) * π[a_i, z_i, N_i] * P[z_i, zp_i]                    
                        T[s_i, sp_i+1] += q * π[a_i, z_i, N_i] * P[z_i, zp_i]                    
        
        # Check that the transition matrix adds up to 1 for each row
        assert(np.allclose(T.sum(1), 1))

        return T

    def find_stationary_distribution(self, T):

        # Find eigenvalue closest to 1
        w, v = sp.sparse.linalg.eigs(T.T, k=1, sigma=1)

        # Find stationary distribution
        λ = np.real(v.flatten())

        return λ/λ.sum()

    def find_aggregate_var(self, g, agrid, zgrid, λ):

        return np.dot(λ, g.reshape(-1, order='F'))

    def find_equilibrium(self, agrid, zgrid, R, w, p):
        find_policies = self.find_policies
        find_transition_matrix = self.find_transition_matrix
        find_stationary_distribution = self.find_stationary_distribution
        find_aggregate_var = self.find_aggregate_var
        N = self.N

        # Find policies
        gc, ga, π = find_policies(agrid, zgrid, R, w, p)

        # Find transition matrix
        T = find_transition_matrix(agrid, zgrid, π, ga)

        # Find stationary distribution
        λ = find_stationary_distribution(T)
        λ_a = np.sum(λ.reshape([agrid.size, zgrid.size], order='F'),1)

        return AiyagariEquilibrium(gc, ga, π, T, λ, λ_a)

class AiyagariEquilibrium:

    def __init__(self, gc, ga, π, T, λ, λ_a):
        self.gc, self.ga, self.π = gc, ga, π
        self.T, self.λ, self.λ_a = T, λ, λ_a


