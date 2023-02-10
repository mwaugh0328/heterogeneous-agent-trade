include("gravity-tools.jl")
include("trade-environment.jl")
using CSV
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


grv_params = gravity_params(L = dflabor.L, dfcntryfix = dfcntryfix)

################################################################
# Run the Gravity regression

trade_cost_type = "waugh"

grvdata = gravity(dftrade, trade_cost_type = trade_cost_type, display = true);

trc = trade_costs(grvdata.dist_coef, grvdata.lang_coef, grvdata.θm)

d = zeros(19,19)

make_trade_costs!(grvdata, d, grv_params, trade_cost_type = trade_cost_type)

T = zeros(19)
W = ones(19)

make_technology!(grvdata, T, W, grv_params)


@time πshares, Φ = eaton_kortum(W, d, T, grv_params.θ)

# Now re-run the gravity regression on the model to 
# see if we recover the proper coeffecients

trademodel = log.(vec(normalize_by_home_trade(πshares, grv_params.Ncntry)'))

dfmodel = DataFrame(trade = trademodel)

filter!(row -> ~(row.trade ≈ 1.0), dfmodel);

filter!(row -> ~(row.trade ≈ 0.0), dfmodel);

dfmodel = hcat(dfmodel, dfcntryfix)

plot(dfmodel.trade, dftrade.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)