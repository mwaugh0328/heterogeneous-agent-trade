include("ha-trade.jl")
using MINPACK
using Plots

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

W_ss = [exp.(sol.x[1]); 1.0]
R_ss = [exp.(sol.x[2]); exp.(sol.x[2])]

# New stuff
tariff = cntry_prm.tariff
σ = 2.0

p = make_p(W_ss[1], TFP, d[1, :], tariff[1, :] )
_, gc_HAT, _, _ = policy_function_itteration(R_ss[1], W_ss[1], p, τ[1], hh_prm)
gc_CES,_,_ = basic_CES(W_ss[1], p, σ, hh_prm)

gc_HAT
gc_CES

plot(gc_HAT[1,:,1],label="HAT")
plot!(gc_CES[1,:,1],label="Basic CES")

plot(gc_HAT[1,:,2],label="HAT")
plot!(gc_CES[1,:,2],label="Basic CES")