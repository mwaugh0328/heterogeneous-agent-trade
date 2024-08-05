include("ha-trade.jl")

using Plots


########################################################################
#
# Parameters
# 


γ = 1.5 # curvatuve on CRRA utility function
σϵ = 0.25 # logit dispersion parameter
Ncntry = 2 # number of countries

# this setups up parameters on the household side
hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92, γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ);

# here are some simple country parametrers
TFP = [1.0; 1.0]

τ = [0.0; 0.0]

L = [1.0; 1.0]

d_ij = 1.745

d = [1.0 d_ij; d_ij 1.0]

# this sets up the country specific paramters
cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP);

# Taking prices as given...
init_R = 1.0
init_W = 1.02
init_p = make_p(init_W, TFP, d[1, :], τ[1, :])

########################################################################
#
# First check
# 

# Check if 'policy_function_itteration' and 'policy_function_itteration_new' give exact same outputs
# still have two versions...
Kga_old, Kgc_old, πprob_old, Tv_old = policy_function_itteration_old(init_R, init_W, init_p, τ[1], hh_prm; tol = 10^-6, Niter = 500)
Kga, Kgc, πprob, Tv = policy_function_itteration(init_R, init_W, init_p, τ[1], hh_prm; tol = 10^-6, Niter = 500)

vec_max(Kga_old, Kga)
vec_max(Kgc_old, Kgc)
vec_max(πprob_old, πprob)
vec_max(Tv_old, Tv)

########################################################################
#
# Second check
# 

# check if we input stationary object (policy function for c and value function), 
# one step iteration gives the same things

# Here we have asset, consumption policies, and choice prob. same, but value function is not
hh_obj = one_step_itteration(Kgc_old, Tv_old, init_R, init_R, init_W, init_p, init_p, τ[1], hh_prm)
vec_max(Kga_old, hh_obj.asset_policy)
vec_max(Kgc_old, hh_obj.cons_policy)
vec_max(πprob_old, hh_obj.πprob)
vec_max(Tv_old, hh_obj.Tv)

#a, b, c = coleman_operator_new(Kgc, Tv, init_R, init_R, init_W, init_p, init_p, τ[1], hh_prm)
#vec_max(b, Tv)

#d,e,f = coleman_operator(Kgc, Tv, init_R, init_W, init_p, τ[1], hh_prm)
#vec_max(e, Tv)