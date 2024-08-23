include("ha-trade.jl")

using MINPACK
using JET


# ##########################################################################
# Parameters
# 

γ = 1.5 # curvatuve on CRRA utility function
σϵ = 0.25 # logit dispersion parameter
Ncntry = 2 # number of countries

hh_prm = household_params();

# here are some simpe country parametrers
TFP = [1.0; 1.0]

τ = [0.0; 0.0]

L = [1.0; 1.0]

d_ij = 1.745

d = [1.0 d_ij; d_ij 1.0]

# this sets up the country specific paramters
cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP);

# ##########################################################################
# Equilibrium
# 

f(x) = world_equillibrium_FG(exp.(x), hh_prm, cntry_prm)
# this world... function is used to construct zero conditions for
# the finacial globalization case

function f!(fvec, x)

    fvec .= f(x)

end


xguess = [1.0; 1.02]

n = length(xguess)
diag_adjust = n - 1

sol = fsolve(f!, log.(xguess), show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-10,
       )

# This solver works very well in general. Spencer Lyon converted old-school minpack
# routines in C into julia

print(sol)

wage = [exp.(sol.x[1]); 1.0]
R = exp.(sol.x[2])

Rsol = [R; R]

Y, tradeflows, A_demand, Gbudget, tradeshare, hh, dist = world_equillibrium(Rsol, wage, τ, hh_prm, cntry_prm, tol_vfi = 1e-10);
# this world_eq...is a core file takes prices and returns a bunch of stuff
# note that hh, dist are objects of dimensiom number of countries, then within
# it has policy functions and distributions state by state

# ##########################################################################
# Transition path
# 

T = 3

TFP = [1.0; 1.0]

τ = [0.0; 0.0]

L = [1.0; 1.0]

tariff = zeros(Ncntry, Ncntry)

d_ij = 1.745
d = [1.0 d_ij; d_ij 1.0]

Rpath = repeat(Rsol, outer = (1,T)) 
Rend = copy(Rsol)

d_path = d.* ones(Ncntry, Ncntry, T+1)

W_path = repeat( [1.0; 1.0], outer = (1,T))
xxx = W_path[:]
Wend = copy(wage)

trp_values = trans_path_values(hh, dist, Rend, Wend, T, τ)

good_market, asset_market = transition_path(xxx, Rpath, d_path, trp_values, hh_prm, cntry_prm)

@report_opt transition_path(xxx, Rpath, d_path, trp_values, hh_prm, cntry_prm)
# this run gives 458 possible errors





