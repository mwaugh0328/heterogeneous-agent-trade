include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("static-trade-environment.jl")
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

df = DataFrame(CSV.File("solution-fg.csv"))

initial_x = [df.wage[2:end]; 1.035]

TFP = df.TFP
L = df.L

Ncntry = size(d)[1]

mdl_prm = world_model_params(Ncntry = Ncntry, Na = 50, 
γ = 1.5, ϕ = 2.0, amax = 5.0, σϵ = 0.25, d = d, TFP = TFP, L = L)

p = (df.wage[1:end] ./ mdl_prm.TFP).*d[1,:]

compute_eq(1.02, df.wage[1], p, mdl_prm)

# trade with interpolations at .13.1
# not multi-threaded on policy function, only value fun
# @time compute_eq(1.02, df.wage[1], p, mdl_prm);
# 76
#   0.776629 seconds (35.50 k allocations: 44.101 MiB)
#   0.005110 seconds (28.68 k allocations: 3.060 MiB)
#   0.782261 seconds (64.23 k allocations: 47.165 MiB)


# this is ha-mc with interpolations lastest version
# # @time hh_end, dist_end = compute_eq(R, W, πrft, p, mdl_prm,
# hh_solution_method = "itteration", stdist_sol_method = "itteration");
# 61
#   0.725138 seconds (183.59 k allocations: 79.649 MiB, 2.20% gc time)




