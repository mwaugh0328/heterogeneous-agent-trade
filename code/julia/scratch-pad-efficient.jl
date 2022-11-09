include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("static-trade-environment.jl")
include("ha-efficient.jl")
using MINPACK
using Plots
using CSV
using DataFrames

####################################################################################
####################################################################################
println(" ")
println(" ")
println("########### computing initial eq ################")
println(" ")

Ncntry = 19

dftrade = DataFrame(CSV.File("ek-trade.csv"))

d = reshape(dftrade.d, Ncntry,Ncntry)

df = DataFrame(CSV.File("solution.csv"))

initial_x = [df.wage[2:end]; df.interest_rate]

#initial_x =[df.TFP[2:end]; 1.02*ones(Ncntry)]

TFP = df.TFP
L = df.L

Ncntry = size(d)[1]

mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, 
γ = 1.5, ϕ = 2.0, amax = 8.0, σϵ = 0.25, d = d, TFP = TFP, L = L)

####################################################################################

initial_x = 0.8.*ones(mdl_prm.Ncntry)

f(x) = efficient_equillibrium(x, mdl_prm)

function f!(fvec, x)

    fvec .= f(x)

end

n = length(initial_x)
diag_adjust = n - 1

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,)

print(sol)

social = compute_efficient(sol.x, mdl_prm)