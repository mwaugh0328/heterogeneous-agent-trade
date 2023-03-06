include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("ha-trade-elasticity.jl")
include("ha-efficient.jl")
include("static-trade-environment.jl")
include("gravity-tools.jl")
using MINPACK
using Plots
using CSV
using DataFrames

####################################################################################
####################################################################################
# This sets up the EK trade data and gravity stuff

dftrade = DataFrame(CSV.File("../../ek-data/ek-data.csv"))

dftrade.trade = parse.(Float64, dftrade.trade)
    # forsome reason, now it reads in as a "String7"
    
dflang = DataFrame(CSV.File("../../ek-data/ek-language.csv"))
    
dflabor = DataFrame(CSV.File("../../ek-data/ek-labor.csv"))
    
filter!(row -> ~(row.trade ≈ 1.0), dftrade);
    
filter!(row -> ~(row.trade ≈ 0.0), dftrade);
    
dftrade = hcat(dftrade, dflang);
    
    #dfcntryfix = select(dftrade,Not("trade"))
dfcntryfix = DataFrame(CSV.File("../../ek-data/ek-cntryfix.csv"))
    # these are the fixed characteristics of each country...


trade_cost_type = "ek"

grvdata = gravity(dftrade, display = true, trade_cost_type = trade_cost_type );

trc = trade_costs(grvdata.dist_coef, grvdata.lang_coef, grvdata.θm)

grv_params = gravity_params(L = dflabor.L, dfcntryfix = dfcntryfix, Ncntry = 19)

####################################################################################
####################################################################################

df = DataFrame(CSV.File("solution-fg.csv"))

L = df.L

Ncntry = size(L)[1]

γ = 1.5
σϵ = 0.25

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

dfparams = DataFrame(CSV.File("current-guess-ek-scale.csv"))
#dfparams = DataFrame(CSV.File("current-guess-log-ek.csv"))

xxx = dfparams.guess

out_moment_vec, Wsol, Rsol, πshare  = calibrate(xxx, grvdata, grv_params, hh_prm, cntry_prm, trade_cost_type = trade_cost_type)

TFP, d = make_country_params(xxx, cntry_prm, grv_params, trade_cost_type = trade_cost_type)

cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP)

####################################################################################
####################################################################################
println(" ")
println(" ")
println("########### computing initial eq ################")
println(" ")

initial_x = 0.65*log.(TFP)

ψ = exp.(0.95*log.(TFP))

f(x) = efficient_equillibrium(exp.(x), ψ , hh_prm, cntry_prm)

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

social = compute_efficient(exp.(sol.x), ψ , hh_prm, cntry_prm)

####################################################################################
####################################################################################

plot(dftrade.trademodel, dftrade.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)



####################################################################################
####################################################################################











# f(x) = world_equillibrium(x, mdl_prm, hh_solution_method = "itteration", stdist_sol_method = "itteration");

# function f!(fvec, x)

#     fvec .= f(x)

# end

# ###################################################################

# n = length(initial_x)
# diag_adjust = n - 1

# sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
#       ml=diag_adjust, mu=diag_adjust,
#       diag=ones(n),
#       mode= 1,
#       tol=1e-5,
#        )

# print(sol)

# Wsol = [1.0; sol.x[1:Ncntry-1]]
# Rsol = sol.x[Ncntry:end]

# Y, tradeflows, A_demand, tradeshare, hh, dist = world_equillibrium(Rsol,
#     Wsol, mdl_prm, hh_solution_method = "itteration");

# # plot(log.(vec(tradeshare)), log.(dftrade.tradesharedata), seriestype = :scatter)

# dftrade_model_data = DataFrame(
#     trademodel = vec(tradeshare),
#     tradedata = dftrade.tradesharedata,
#     trade_efficient = vec(social.tradeshare),
#     tradecost = dftrade.d,
#     importer_index = dftrade.importer_index,
#     exporter_index = dftrade.exporter_index
#      );

# # CSV.write("../../notebooks/trade_model_data.csv", dftrade_model_data)

# # hh_df = make_hh_dataframe(dist, hh, 19, Rsol, Wsol, mdl_prm)

# # CSV.write("../../notebooks/household_data_pre.csv", hh_df)



