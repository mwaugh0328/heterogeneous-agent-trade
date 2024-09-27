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
# Find initial Equilibrium

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

Wint = [exp.(sol.x[1]); 1.0]

Rint = [exp.(sol.x[2]); exp.(sol.x[2])]

Y, tradeflows, A_demand, Gbudget, tradeshare, hh, dist_int = world_equillibrium(Rint, Wint, τ, hh_prm, cntry_prm, tol_vfi = 1e-10);
# this world_eq...is a core file takes prices and returns a bunch of stuff
# note that hh, dist are objects of dimensiom number of countries, then within
# it has policy functions and distributions state by state

# ##########################################################################
# Find ending Equilibrium

d_ij_end = 1.65

d_end = [1.0 d_ij_end; d_ij 1.0]

cntry_prm = country_params(cntry_prm, d = d_end);

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

Wend = [exp.(sol.x[1]); 1.0]

Rsol_end = [exp.(sol.x[2]); exp.(sol.x[2])]

Y, tradeflows, A_demand, Gbudget, tradeshare, hh_end, dist_end = world_equillibrium(Rsol_end, Wend, τ, hh_prm, cntry_prm, tol_vfi = 1e-10);
# this world_eq...is a core file takes prices and returns a bunch of stuff
# note that hh, dist are objects of dimensiom number of countries, then within
# it has policy functions and distributions state by state

# ##########################################################################
# Transition path
# 

T = 150

Rpath_check = repeat(Rint, outer = (1,T-1)) # in PI this is length T?
Rend_check = copy(Rint)

d_path_check = d .* ones(Ncntry, Ncntry, T+1)

W_path_check = repeat( Wint, outer = (1,T))
xxx_check = W_path_check[:]
Wend_check = copy(Wint)

trp_values_check = trans_path_values(hh, dist_int, Rint, Rend_check, Wend_check, T, τ)

good_market_check, asset_market_check = transition_path(xxx_check, Rpath_check, d_path_check, trp_values_check, hh_prm, cntry_prm)

good_market_check
asset_market_check

################################################

Rpath = repeat(Rsol_end, outer = (1,T-1)) # in PI this is length T?
Rend = copy(Rsol_end)

d_path = zeros(Ncntry, Ncntry, T+1)

d_path[:,:,1] = [1.0 1.70; d_ij 1.0]
d_path[:,:,2] = [1.0 1.68; d_ij 1.0]
d_path[:,:,3:end] = d_end .* ones(Ncntry, Ncntry, T-1)

W_path = repeat( Wend, outer = (1,T))
xxx = W_path[:]
Wend = copy(Wend)

trp_values = trans_path_values(hh_end, dist_int, Rint, Rend, Wend, T, τ)

good_market, asset_market = transition_path(xxx, Rpath, d_path, trp_values, hh_prm, cntry_prm)


# transition_path_only_assetmarket(Rpath, W_path[:], d_path, trp_values, hh_prm, cntry_prm)

# @report_opt transition_path(xxx, Rpath, d_path, trp_values, hh_prm, cntry_prm)
# this run gives 458 possible errors


# TFP = [1.0; 1.0]

# τ = [0.0; 0.0]

# L = [1.0; 1.0]

# tariff = zeros(Ncntry, Ncntry)





