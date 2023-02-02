include("gravity-tools.jl")
include("trade-environment.jl")
using CSV
using DataFrames
using Plots
using MINPACK

################################################################
# builds the EK dataset

dftrade = DataFrame(CSV.File("../../ek-data/ek-data.csv"))

dflang = DataFrame(CSV.File("../../ek-data/ek-language.csv"))

dflabor = DataFrame(CSV.File("../../ek-data/ek-labor.csv"))

filter!(row -> ~(row.trade ≈ 1.0), dftrade);

filter!(row -> ~(row.trade ≈ 0.0), dftrade);

dftrade = hcat(dftrade, dflang);

#dfcntryfix = select(dftrade,Not("trade"))
dfcntryfix = DataFrame(CSV.File("../../ek-data/ek-cntryfix.csv"))
# these are the fixed characteristics of each country...
# not trade flows

################################################################
# Run the Gravity regression

grvdata = gravity(dftrade, display = true);


# ################################################################
# # Recover the trade costs and technology parameters

θ = 4.0

d = zeros(19,19)
T = zeros(19)
W = ones(19)

make_trade_costs!(dfcntryfix, grvdata, d, θ)

make_technology!(grvdata, T, W, θ)

#@time gravity!(dftrade, d, T, W, θ)

################################################################
# Feed into the model and then compare with data

@time πshares, Φ = eaton_kortum(W, d, T, θ)

# Now re-run the gravity regression on the model to 
# see if we recover the proper coeffecients

trademodel = log.(vec(normalize_by_home_trade(πshares)'))

dfmodel = DataFrame(trade = trademodel)

filter!(row -> ~(row.trade ≈ 1.0), dfmodel);

filter!(row -> ~(row.trade ≈ 0.0), dfmodel);

dfmodel = hcat(dfmodel, dfcntryfix)

grv = gravity(dfmodel, display = true);

plot(dfmodel.trade, dftrade.trade, seriestype = :scatter)

################################################################
# Feed into the model and then compare with data


