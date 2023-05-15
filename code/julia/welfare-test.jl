include("ha-trade.jl")
include("ha-trade-welfare.jl")

using MINPACK
using Plots
using CSV
using DataFrames


γ = 1.05 # curvatuve on CRRA utility function
σϵ = 0.25 # logit dispersion parameter
Ncntry = 2 # number of countries

# this setups up parameters on the household side
hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92, γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ);


# here are some simpe country parametrers
TFP = [1.0; 1.0]

τ = [0.0; 0.0]

L = [1.0; 1.0]

d_ij = 2.05

d = [1.0 d_ij; d_ij 1.0]

# this sets up the country specific paramters
cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP);

#####################################################################################
#####################################################################################
# now compute initial eq.

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

Y, tradeflows, A_demand, Gbudget, tradeshare, hh, dist = world_equillibrium(Rsol, wage, τ, hh_prm, 
                                cntry_prm, tol_vfi = 1e-10);

#####################################################################################
#####################################################################################
# now compute counterfact eq.

Δd = 0.10

d = [1.0 d_ij * (1 - Δd)  ; d_ij *(1 - Δd)  1.0]

# this sets up the country specific paramters
Δd_cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP);

#####################################################################################
#####################################################################################
# now new eq.

f(x) = world_equillibrium_FG(exp.(x), hh_prm, Δd_cntry_prm)
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

Δwage = [exp.(sol.x[1]); 1.0]
ΔR = exp.(sol.x[2])

ΔRsol = [ΔR; ΔR]


ΔY, Δtradeflows, ΔA_demand, ΔGbudget, Δtradeshare, Δhh, Δdist = world_equillibrium(ΔRsol, Δwage, τ, hh_prm, 
                                                Δd_cntry_prm, tol_vfi = 1e-10);

#####################################################################################
#####################################################################################

cntry = 1

p = make_p(wage, cntry_prm.TFP, cntry_prm.d[cntry, :], cntry_prm.tariff[cntry, :] )

τeqv =  eq_variation(Rsol[cntry], wage[cntry], p, Δhh[cntry], dist[cntry].state_index, hh_prm)

pcons, xpcons = bilateral_consumption(Rsol[cntry], wage[cntry], hh, cntry, hh_prm)

∂W, ∂logW = welfare_by_state(hh, Δhh, cntry, hh_prm.σϵ)