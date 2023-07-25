import numpy as np
import scipy as sp
from copy import copy
from hat_support import *

class HATmodel:

    def __init__(self, β, γ, φ, σ_ε, N, P, L, TFP, τ, d, agrid, zgrid, ψ):
        self.β, self.γ, self.φ = β, γ, φ
        self.σ_ε, self.N, self.P = σ_ε, N, P
        self.L, self.TFP, self.τ, self.d = L, TFP, τ, d 
        self.agrid, self.zgrid, self.ψ = agrid, zgrid, ψ

    def utility(self, c):
        γ = self.γ
        return ((1-γ)**(-1))*(c**(1-γ))

    def u_c(self, c):
        γ = self.γ        
        return c**(-γ)

    def u_c_inv(self, x):
        γ = self.γ        
        return x**(-1/γ) 

    def find_p(self, w):
        TFP, d, τ = self.TFP, self.d, self.τ
        return (w / TFP) * d * (1 + τ)

    def find_Y(self):
        TFP, L = self.TFP, self.L
        return TFP*L

    def change_d(self, d):
        self.d = d

    def coleman(self, x, R, w):
        u_c, u_c_inv = self.u_c, self.u_c_inv
        β, φ, P, N = self.β, self.φ, self.P, self.N
        agrid, zgrid = self.agrid, self.zgrid
        find_value_function = self.find_value_function
        find_probabilities = self.find_probabilities
        find_p = self.find_p

        # Find prices
        p = find_p(w)

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

        v1 = find_value_function(v, gc1, ga1)
        x1 = np.concatenate([gc1, v1])

        return x1 

    def find_value_function(self, v, gc, ga):
        find_log_sum, utility = self.find_log_sum, self.utility
        P, N = self.P, self.N
        agrid, zgrid = self.agrid, self.zgrid
        β, ψ = self.β, self.ψ

        for N_i in range(N):
            for z_i, z in enumerate(zgrid):
                for a_i, a in enumerate(agrid): 

                    # Compute Expectations 
                    Ev = 0
                    if ga[a_i, z_i, N_i] <= agrid[0]:
                        Ev += P[z_i, :] @ find_log_sum(v[0, :, :], ψ[0, :, :], zgrid)
                    elif ga[a_i, z_i, N_i] >= agrid[-1]:
                        Ev += P[z_i, :] @ find_log_sum(v[agrid.shape[0]-1, :, :], ψ[agrid.shape[0]-1, :, :], zgrid)
                    else: 
                        
                        idx = np.where(ga[a_i, z_i, N_i]>=agrid)[0][-1]
                        q = ((ga[a_i, z_i, N_i] - agrid[idx]) / (agrid[idx+1] - agrid[idx]))
                        Ev += (1-q) * P[z_i, :] @ find_log_sum(v[idx, :, :], ψ[idx, :, :], zgrid)
                        Ev += q * P[z_i, :] @ find_log_sum(v[idx+1, :, :], ψ[idx+1, :, :], zgrid)

                    v[a_i, z_i, N_i] = utility(gc[a_i, z_i, N_i]) + β*Ev

        return v

    def find_policies(self, R, w):
        u_c, u_c_inv = self.u_c, self.u_c_inv
        coleman, find_probabilities = self.coleman, self.find_probabilities      
        φ, β, N = self.φ, self.β, self.N
        agrid, zgrid = self.agrid, self.zgrid
        fixed_point_iteration = self.fixed_point_iteration
        find_p = self.find_p

        # Find prices
        p = find_p(w)

        # Initial guess for consumption policy
        gc0_vec = np.linspace(0.1, agrid[-1], agrid.shape[0])
        gc0 = np.repeat(np.repeat(gc0_vec[:,np.newaxis], zgrid.shape[0], axis=1)[:,:,np.newaxis], N, axis=2) 
        v0 = (-1)*np.ones((agrid.shape[0], zgrid.shape[0], N))/(1-β)
        x0 = np.concatenate([gc0, v0])

        # Find fixed point of coleman operator
        fun = lambda x: coleman(x, R, w)     
        x = fixed_point_iteration(fun, x0, agrid)

        gc = x[0:agrid.shape[0], :, :]
        v = x[agrid.shape[0]:, :, :]

        # Find labor and asset policies
        cih = zgrid*w + R*agrid[:,np.newaxis]
        ga = np.maximum(cih[:,:,np.newaxis] - gc*p, φ)       

        # Find probabilities
        π = find_probabilities(v)

        return gc, ga, π, v

    def fixed_point_iteration(self, fun, x0, agrid, max_iter = 1000, error=1e-8):

        x = x0
        gc = x[0:agrid.shape[0], :, :]
        for i in range(max_iter):
            x1 = fun(x)
            gc1 = x1[0:agrid.shape[0], :, :]

            if np.max(np.abs(gc1-gc)) < error:
                return x1
            else:
                x = x1
                gc = gc1

        print('Fixed point iteration did not converse!')

    def find_probabilities(self, v):
        σ_ε = self.σ_ε
        ψ = self.ψ

        v_temp = ((v + ψ).T - np.max(v + ψ, axis=2).T).T
        return (np.exp(v_temp/σ_ε).T/np.exp(v_temp/σ_ε).sum(2).T).T

    def find_log_sum(self, vj, ψ_j, zgrid):
        σ_ε = self.σ_ε

        vj_max = np.max(vj + ψ_j, axis=1)
        log_sum = np.array([σ_ε*np.log(np.sum(np.exp((vj[z_i, :] + ψ_j[z_i, :] - vj_max[z_i])/σ_ε)))+vj_max[z_i] for z_i in range(len(zgrid))])

        return log_sum

    def find_transition_matrix(self, π, ga):
        """
        Transition matrix for the vector s = [(a1,z1),...,(aM,z1),(a1,z2),...,(aM,z2),...(aM,zN)]'        
        """
        P, N = self.P, self.N
        agrid, zgrid = self.agrid, self.zgrid

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

    def solve_model(self, R, w, cntry_idx=0):
        find_policies = self.find_policies
        find_transition_matrix = self.find_transition_matrix
        find_stationary_distribution = self.find_stationary_distribution
        find_p, find_aggregates = self.find_p, self.find_aggregates

        # Find places
        p = find_p(w)

        # Find policies
        gc, ga, π, v = find_policies(R, w)

        # Find transition matrix
        T = find_transition_matrix(π, ga)

        # Find stationary distribution
        λ = find_stationary_distribution(T)

        # Find aggregate quantities
        A_i, Y_i, C_i, C_ij, G_i = find_aggregates(gc, ga, π, p, λ, cntry_idx)

        return A_i, Y_i, C_i, C_ij, G_i, gc, ga, π, v, T, λ 

    def find_aggregates(self, gc, ga, π, p, λ, cntry_idx):
        agrid, zgrid = self.agrid, self.zgrid
        TFP, L, d, τ = self.TFP, self.L, self.d, self.τ
        N = self.N

        # Compute aggregate assets
        λ_a = np.sum(λ.reshape([agrid.size, zgrid.size], order='F'), 1)
        A = λ_a @ agrid 

        # Output
        λ_e = np.sum(λ.reshape([agrid.size, zgrid.size], order='F'), 0)
        Y = p[cntry_idx] * TFP * L * λ_e @ zgrid

        # Total Consumption and consumption by country
        C = L * np.reshape(np.sum(p * gc * π, 2), -1, order='F') @ λ 
        C_j = L *  np.reshape(p * gc * π, (agrid.shape[0]*zgrid.shape[0], N), order='F').T @ λ

        # Government
        G = L * τ 

        return A, Y, C, C_j, G


class MulticountryHATmodel:

    def __init__(self, β, γ, φ, σ_ε, N, P, L, TFP, τ, d, agrid, zgrid, ψ_slope, numeraire_cntry = 0):
        self.β, self.γ, self.φ = β, γ, φ
        self.σ_ε, self.N, self.P = σ_ε, N, P
        self.L, self.TFP, self.τ, self.d = L, TFP, τ, d
        self.τ, self.d = τ, d
        self.agrid, self.zgrid = agrid, zgrid
        self.ψ_slope = ψ_slope
        self.numeraire_cntry = numeraire_cntry 

    def solve_model_FG(self, R, w):
        β, γ, φ, σ_ε = self.β, self.γ, self.φ, self.σ_ε
        N, P, L, TFP = self.N, self.P, self.L, self.TFP
        τ, d = self.τ, self.d
        agrid, zgrid = self.agrid, self.zgrid
        find_trade_stats = self.find_trade_stats 
        find_ψ = self.find_ψ

        A_i = np.empty(N)
        Y_i = np.empty(N)
        C_i = np.empty(N)
        G_i = np.empty(N)
        C_ij = np.empty((N, N))
        gc_i = np.empty((agrid.shape[0], zgrid.shape[0], N, N)) 
        ga_i = np.empty((agrid.shape[0], zgrid.shape[0], N, N)) 
        v_i = np.empty((agrid.shape[0], zgrid.shape[0], N, N)) 
        π_i = np.empty((agrid.shape[0], zgrid.shape[0], N, N)) 
        T_i = np.empty((agrid.shape[0]*zgrid.shape[0],agrid.shape[0]*zgrid.shape[0], N)) 
        λ_i = np.empty((agrid.shape[0]*zgrid.shape[0], N)) 
        for i in range(N):
            ψ = find_ψ(i)
            m_i = HATmodel(β, γ, φ, σ_ε, N, P, L[i], TFP[i], τ[i], d[i,:], agrid, zgrid, ψ)    
            A_i[i], Y_i[i], C_i[i], C_ij[i, :], G_i[i], gc_i[:,:,:,i], ga_i[:,:,:,i], π_i[:,:,:,i], v_i[:,:,:,i], T_i[:,:,i], λ_i[:,i] = m_i.solve_model(R, w[i], i)

        # Find total exports and imports
        M_i, X_i, tradeshares = find_trade_stats(C_ij)

        return A_i, Y_i, C_i, C_ij, G_i, M_i, X_i, tradeshares, gc_i, ga_i, π_i, v_i, T_i, λ_i

    def objective_FG(self, x):
        numeraire_cntry = self.numeraire_cntry
        solve_model_FG = self.solve_model_FG
        N = self.N

        x = np.exp(x)
        R = x[0]
        w = np.zeros(N)
        mask = np.zeros(N, dtype=bool)
        mask[numeraire_cntry] = True
        w[mask] = 1
        w[mask==False] = x[1:]       

        A_i, Y_i, C_i, C_ij, G_i, M_i, X_i, tradeshare, gc_i, ga_i, π_i, v_i, T_i, λ_i = solve_model_FG(R, w)
        Y_i1 = C_i + G_i + X_i - M_i 

        resid = np.zeros(N)
        resid[0] = A_i.sum()
        resid[1:] = Y_i[mask==False] - Y_i1[mask==False] 

        return resid

    def find_trade_stats(self, C_ij):
        N = self.N

        X_i = np.empty(N)
        M_i = np.empty(N)
        tradeshares = np.empty((N, N))
        for i in range(N):
            M_i[i] = C_ij[i, :].sum() - C_ij[i, i]
            X_i[i] = C_ij[:, i].sum() - C_ij[i, i]        
            tradeshares[i,:] = C_ij[i, :]/C_ij[i, :].sum() 

        return M_i, X_i, tradeshares

    def find_equilibrium_FG(self, x0):
        objective_FG = self.objective_FG 
        solve_model_FG = self.solve_model_FG
        numeraire_cntry = self.numeraire_cntry
        N = self.N  

        out = sp.optimize.root(objective_FG, np.log(x0))
        out = np.exp(out.x)

        R = out[0]
        w = np.zeros(N)
        mask = np.zeros(N, dtype=bool)
        mask[numeraire_cntry] = True
        w[mask] = 1
        w[mask==False] = out[1:]    

        A_i, Y_i, C_i, C_ij, G_i, M_i, X_i, tradeshares, gc_i, ga_i, π_i, v_i, T_i, λ_i = solve_model_FG(R, w) 

        return Equilibrium(R, w, A_i, Y_i, C_i, C_ij, G_i, M_i, X_i, tradeshares, gc_i, ga_i, π_i, v_i, T_i, λ_i) 

    def find_log_c(self, m_i, R, w, log_d):

        # Change d
        m_i1 = copy(m_i)
        m_i1.change_d(np.exp(log_d))

        # Find policies
        gc, ga, π, v = m_i1.find_policies(R, w)

        return np.log(gc), np.log(π) 

    def find_trade_elasticities(self, R, w, i, h=1e-8):
        β, γ, φ, σ_ε = self.β, self.γ, self.φ, self.σ_ε
        N, P, L, TFP = self.N, self.P, self.L, self.TFP
        τ, d, N = self.τ, self.d, self.N
        agrid, zgrid = self.agrid, self.zgrid
        find_log_c = self.find_log_c

        # Model instance for country i
        m_i = HATmodel(β, γ, φ, σ_ε, N, P, L[i], TFP[i], τ[i], d[i,:], agrid, zgrid)    

        # Function we want to take the derivative of
        fn = lambda log_d: find_log_c(m_i, R, w[i], log_d)

        θc = np.zeros((agrid.shape[0], zgrid.shape[0], N))
        θπ = np.zeros((agrid.shape[0], zgrid.shape[0], N))
        for j in range(N):
            hvec = np.zeros(N)
            hvec[j] = h
            
            gc1, π1 = fn(np.log(d[i,:]) + hvec) 
            gc2, π2 = fn(np.log(d[i,:])) 

            θc[:,:,j] = h**(-1)*(gc1[:,:,j] - gc2[:,:,j])
            θπ[:,:,j] = h**(-1)*(π1[:,:,j] - π2[:,:,j]) 

        return θc, θπ 

    def find_ψ(self, cntry):
        agrid, zgrid, N = self.agrid, self.zgrid, self.N
        ψ_slope = self.ψ_slope 

        slope = ψ_slope * np.log(zgrid)

        ψ = np.zeros((agrid.shape[0], zgrid.shape[0], N))
        for a_i, a in enumerate(agrid):
            for z_i, z in enumerate(zgrid):
                for N_i in range(N):
                    if cntry == N_i:
                        ψ[a_i, z_i, N_i] = 0. + slope[z_i] 

        return ψ 

            
class Equilibrium:

    def __init__(self, R, w, A_i, Y_i, C_i, C_ij, G_i, M_i, X_i, tradeshares, gc_i, ga_i, π_i, v_i, T_i, λ_i):
        self.R, self.w = R, w
        self.A_i, self.Y_i, self.C_i = A_i, Y_i, C_i 
        self.C_ij, self.G_i, self.M_i = C_ij, G_i, M_i 
        self.X_i, self.tradeshares = X_i, tradeshares 
        self.gc_i, self.ga_i, self.π_i, self.v_i, self.T_i, self.λ_i = gc_i, ga_i, π_i, v_i, T_i, λ_i


class CEV:

    def __init__(self, m1, m2, eqm1, eqm2):
        self.m1, self.m2, self.eqm1, self.eqm2 = m1, m2, eqm1, eqm2 
        self.N = m1.N
        self.p1, self.p2 = m1.find_p(eqm1.w), m2.find_p(eqm2.w)

    def objective(self, CEV, i):
        m1, m2, eqm1, eqm2  = self.m1, self.m2, self.eqm1, self.eqm2 
        p1, p2 = self.p1, self.p2
        N = self.N

