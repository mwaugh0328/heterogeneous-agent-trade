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

df = DataFrame(CSV.File("solution.csv"))

initial_x = [df.wage[2:end]; df.interest_rate]

#initial_x =[df.TFP[2:end]; 1.02*ones(Ncntry)]

TFP = df.TFP
L = df.L

γ = 1.5

Ncntry = size(d)[1]

mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, 
γ = γ, ϕ = 2.0, amax = 8.0, σϵ = 0.25, d = d, TFP = TFP, L = L)

#####################################################################################
println(" ")
println(" ")
println("########### computing trade flows, fixed prices ################")
println(" ")

Wsol = df.wage
Rsol = df.interest_rate

Y, tradeflows, A_demand, tradeshare, hh, dist = world_equillibrium(Rsol,
    Wsol, mdl_prm, hh_solution_method = "itteration");

plot(log.(vec(tradeshare)), log.(dftrade.tradesharedata), seriestype = :scatter)

#####################################################################################
println(" ")
println(" ")
println("########### computing Δ trade flows, fixed prices ################")
println(" ")

Δ_d = 0.01

d_prime =  1.0 .+ (d .- 1.0).*(1.0 - Δ_d)

Δ_mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, 
γ = γ, ϕ = 2.0, amax = 8.0, σϵ = 0.25, d = d_prime, TFP = TFP, L = L)

Δ_Y, Δ_tradeflows, Δ_A_demand, Δ_tradeshare, Δ_hh, Δ_dist = world_equillibrium(Rsol,
Wsol, Δ_mdl_prm, hh_solution_method = "itteration");

global_trade_elasticity =  (log.(Δ_tradeshare ./ diag(Δ_tradeshare)) .- log.(tradeshare ./ diag(tradeshare))) ./ (log.(d_prime) .- log.(d))


df





dftrade_model_data = DataFrame(
    importer_index = dftrade.importer_index,
    exporter_index = dftrade.exporter_index,
    elasticity = vec(global_trade_elasticity),
    trademodel = log.(vec(tradeshare)),
    tradedata = log.(dftrade.tradesharedata)
     );

CSV.write("trade_model_data.csv", dftrade_model_data)