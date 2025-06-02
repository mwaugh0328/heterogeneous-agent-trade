# This file computes HAT model with log utility with tariff to compare results with Armington model
include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
using MINPACK

#########################################################################
#
# Set Parameters
#

γ = 1.0 # log utility
σϵ = 0.25 # logit dispersion parameter
Ncntry = 2 # number of countries

R = [1.014, 1.014] # interest rate
W = [1.0; 1.0]
p = [0.9; 1.1]
τ = [0.0; 0.0]

# here are some simpe country parametrers
TFP = [1.0; 1.0]

τ = [0.01; 0.01]

L = [1.0; 1.0]

d_ij = 1.65

d = [1.0 d_ij; d_ij 1.0]    

tariff = [0.0 0.25; 0.25 0.0]

# this sets up the country specific paramters
cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP, tariff = tariff)

#########################################################################
#
# Solve household problem and aggregation (HAT with log utility)
#

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92, γ = γ, ϕ = 0.5, amax = 3.0, σϵ = 0.25);

Kga, Kgc, πprob, Tv = policy_function_itteration(R[1], W[1], p, τ[1], hh_prm; tol = 10^-10, Niter = 1000)

hh1, dist1 = compute_eq(R[1], W[1], p, τ[1], hh_prm)

hh, dist, output, tradestats = world_equillibrium(R, W, τ, hh_prm, cntry_prm )

#########################################################################

f(x) = world_equillibrium_FG_tariff((x), hh_prm, cntry_prm)
# this world... function is used to construct zero conditions for
# the finacial globalization case

function f!(fvec, x)

    fvec .= f(x)

end

xguess = [1.0; 0.0; 0.0; 1.02]

n = length(xguess)
diag_adjust = n - 1

sol = fsolve(f!, (xguess), show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-10,
       )


Y, tradeflows, tradeflows_net_tariff, A_demand, Gbudget, tradeshare, hh, dist = world_equillibrium_FG_tariff(sol.x, hh_prm, cntry_prm, display = true)
tradeflows
tradeflows_net_tariff
tradeshare