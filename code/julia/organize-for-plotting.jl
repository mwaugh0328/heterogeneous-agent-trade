include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
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

initial_prices = [df.wage[1:end-1]; 1.00]

L = df.L

Ncntry = size(L)[1]

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = 1.5, ϕ = 0.5, amax = 8.0, σϵ = 0.25)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

dfparams = DataFrame(CSV.File("current-guess-ek-scale.csv"))

xxx = dfparams.guess

out_moment_vec, Wsol, Rsol, πshare  = calibrate(xxx, grvdata, grv_params, hh_prm, cntry_prm, trade_cost_type = trade_cost_type)

####################################################################################
####################################################################################

trademodel = log.(vec(normalize_by_home_trade(πshare, grv_params.Ncntry)'))

dfmodel = DataFrame(trademodel = trademodel)

filter!(row -> ~(row.trademodel ≈ 1.0), dfmodel);

filter!(row -> ~(row.trademodel ≈ 0.0), dfmodel);

dftrade = hcat(dftrade, dfmodel);

CSV.write("model-data-trade.csv", dfguess)

plot(dftrade.trademodel, dftrade.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)
     