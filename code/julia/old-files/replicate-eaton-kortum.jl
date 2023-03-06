include("gravity-tools.jl")
include("trade-environment.jl")
using CSV
using DataFrames
using Plots
using MINPACK

################################################################
# builds the EK dataset

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
# not trade flows


grv_params = gravity_params(L = dflabor.L, dfcntryfix = dfcntryfix)

################################################################
# Run the Gravity regression

grvdata = gravity(dftrade, display = true);


# ################################################################
# # Recover the trade costs and technology parameters

θ = 4.0

d = zeros(19,19)
T = zeros(19)
W = ones(19)

make_trade_costs!(grvdata, d, grv_params)

make_technology!(grvdata, T, W, grv_params)

################################################################
# Feed into the model and then compare with data


@time πshares, Φ = eaton_kortum(W, d, T, θ)

# Now re-run the gravity regression on the model to 
# see if we recover the proper coeffecients

trademodel = log.(vec(normalize_by_home_trade(πshares, grv_params.Ncntry)'))

dfmodel = DataFrame(trade = trademodel)

filter!(row -> ~(row.trade ≈ 1.0), dfmodel);

filter!(row -> ~(row.trade ≈ 0.0), dfmodel);

dfmodel = hcat(dfmodel, dfcntryfix)

grv = gravity(dfmodel, display = true);


plot(dfmodel.trade, dftrade.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)

################################################################
# Feed into the model and then compare with data


