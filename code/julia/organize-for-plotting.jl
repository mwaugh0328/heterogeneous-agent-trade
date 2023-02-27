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

γ = 1.5
σϵ = 0.25

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

#dfparams = DataFrame(CSV.File("current-guess-ek-scale.csv"))
dfparams = DataFrame(CSV.File("current-guess-log-ek.csv"))

xxx = dfparams.guess

out_moment_vec, Wsol, Rsol, πshare  = calibrate(xxx, grvdata, grv_params, hh_prm, cntry_prm, trade_cost_type = trade_cost_type)

TFP, d = make_country_params(xxx, cntry_prm, grv_params, trade_cost_type = trade_cost_type)

####################################################################################
####################################################################################

trademodel = log.(vec(normalize_by_home_trade(πshare, grv_params.Ncntry)'))

dfmodel = DataFrame(trademodel = trademodel)

filter!(row -> ~(row.trademodel ≈ 1.0), dfmodel);

filter!(row -> ~(row.trademodel ≈ 0.0), dfmodel);

dftrade = hcat(dftrade, dfmodel);

#CSV.write("model-data-trade.csv", dftrade)
CSV.write("log-model-data-trade.csv", dftrade)

plot(dftrade.trademodel, dftrade.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)

####################################################################################
####################################################################################
# this sets things up to look at elasticities


cntry = 19

agrid = make_agrid(hh_prm, TFP[cntry])

foo = household_params(hh_prm, agrid = agrid, TFP = TFP[cntry], σϵ = σϵ*(TFP[cntry]^(1.0 - γ)))

#the way this grid is setup seems to work

p = (Wsol[1:end] ./ TFP).*d[cntry,:]

hh = solve_household_problem(Rsol[cntry], Wsol[cntry], p, foo)

dist = make_stationary_distribution(hh, foo)

ϵπ = similar(hh.πprob)

ϵc = similar(hh.πprob)

@time make_ϵ!(ϵπ, ϵc, hh.cons_policy, hh.Tv, Rsol[cntry], Wsol[cntry], p, foo; points = 3, order = 1);