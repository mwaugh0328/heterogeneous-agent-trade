include("ha-trade.jl")
using MINPACK
using JET
using Plots


# ##########################################################################
# Parameters
# 

γ = 1.5 # curvatuve on CRRA utility function
σϵ = 0.25 # logit dispersion parameter
Ncntry = 2 # number of countries

hh_prm = household_params(γ = γ, σϵ = σϵ, Ncntry= Ncntry, Na = 200);

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
Wint = Wint ./ ( sum(Wint, dims = 1) / Ncntry ) # Need to be consistent with numeriare 

Rint = [exp.(sol.x[2]); exp.(sol.x[2])]

Y, tradeflows, tradeflows_net_tariff, A_demand, Gbudget, tradeshare_int, hh, dist_int = world_equillibrium(Rint, Wint, τ, hh_prm, cntry_prm, tol_vfi = 1e-10);
# this world_eq...is a core file takes prices and returns a bunch of stuff
# note that hh, dist are objects of dimensiom number of countries, then within
# it has policy functions and distributions state by state

# ##########################################################################
##########################################################################
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
Wend = Wend ./ ( sum(Wend, dims = 1) / Ncntry ) # Need to be consistent with numeriare 

Rend = [exp.(sol.x[2]); exp.(sol.x[2])]

Y, tradeflows, tradeflow_net_tariff, A_demand, Gbudget, tradeshare, hh_end, dist_end = world_equillibrium(Rend, Wend, τ, hh_prm, cntry_prm, tol_vfi = 1e-10);
# this world_eq...is a core file takes prices and returns a bunch of stuff
# note that hh, dist are objects of dimensiom number of countries, then within
# it has policy functions and distributions state by state

# ##########################################################################
##########################################################################
# Transition path

T = 50

Rpath = 1.015*ones(T-1)

d_path = zeros(Ncntry, Ncntry, T+1)

# d_path[:,:,1] = [1.0 1.72; d_ij 1.0]
# d_path[:,:,2] = [1.0 1.71; d_ij 1.0]

d_path = d_end .* ones(Ncntry, Ncntry, T+1)

W_path = 1.0.*ones(T)

xxx = vcat(Rpath, W_path)

trp_values = trans_path_values(hh_end, dist_int, Rint, Rend, Wend, T, τ)

out = transition_path(xxx, d_path, trp_values, hh_prm, cntry_prm)

# @time good_market, asset_market = transition_path(xxx, Rpath, d_path, trp_values, hh_prm, cntry_prm)

##########################################################################
##########################################################################

g(x) = transition_path(( x ), d_path, trp_values, hh_prm, cntry_prm)

function g!(fvec, x)

    fvec .= g(x)

end

initial_x = (xxx)

n = length(initial_x)
diag_adjust = n - 1

sol = fsolve(g!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

print(sol);  

godos_market, asset_market, hhpath = transition_path(( sol.x ), d_path, trp_values, hh_prm, cntry_prm, display = true)

vint, vnew, evτ, income = one_time_asset(hh[1], hhpath[1,1], Rint[1], Wint[1], 1.0, hh_prm)



