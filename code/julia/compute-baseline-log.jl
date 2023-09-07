include("ha-trade.jl")
using MINPACK
using Plots
using CSV
using DataFrames
using StatsBase

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
# Compute the EQ at the gravity parameters

# dfparams = DataFrame(CSV.File("current-guess-145-30-725.csv"))
dfparams = DataFrame(CSV.File("./calibration-files/current-guess-log-24.csv"))
xxx = dfparams.guess[1:end]

L = dflabor.L

Ncntry = size(L)[1]

γ = 1.0
σϵ = 0.2369
ψslope = 0.0


hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ, ψslope = ψslope)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

# dfparams = DataFrame(CSV.File("./calibration-files/current-guess-ek-quality60.csv"))

R = 1.01

out_moment_vec, Wsol, β, πshare  = calibrate(xxx, R, grvdata, grv_params, hh_prm, 
                            cntry_prm, trade_cost_type = trade_cost_type)[1:4]

TFP, d = make_country_params(xxx, cntry_prm, grv_params, trade_cost_type = trade_cost_type)

cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP)

hh_prm = household_params(hh_prm, β = β)

Rsol = ones(Ncntry)*R

dfWsol = DataFrame(guess = [Wsol[1:18] ; β])   

# CSV.write("./calibration-files/current-wage-guess-log-22.csv", dfWsol)

####################################################################################
####################################################################################

Y, tradeflows, A_demand, Gbudget, tradeshare, hh, dist = world_equillibrium(Rsol, Wsol, hh_prm, cntry_prm, tol_vfi = 1e-10);

# This is a Plot test to make sure this is doing what I think it is

trademodel = log.(vec(normalize_by_home_trade(tradeshare, grv_params.Ncntry)'))

dfplot = DataFrame(trademodel = trademodel)

filter!(row -> ~(row.trademodel ≈ 1.0), dfplot);

filter!(row -> ~(row.trademodel ≈ 0.0), dfplot);

dfplot = hcat(dftrade, dfplot);

plot(dfplot.trademodel, dfplot.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false,
    show = true)

rootfile = "../../notebooks/output/"

root = rootfile*"log-model-data-trade.csv"

CSV.write(root, dfplot)

# ####################################################################################
# ####################################################################################
# # Let's construct bilateral trade elasticities

cntry = 19 # this is the country I'll look at

p = make_p(Wsol[1:end], TFP, d[cntry, :], cntry_prm.tariff[cntry, :] )
# prices from the perspective of those in that country

agrid = make_agrid(hh_prm, TFP[cntry])

ψ = make_ψ(cntry, ψslope.*TFP[cntry].^(1.0 - γ), hh_prm)

agrid = make_agrid(hh_prm, TFP[cntry])

foo_hh_prm = household_params(hh_prm, agrid = agrid, 
TFP = TFP[cntry], L = L[cntry], σϵ = σϵ*(TFP[cntry]^(1.0 - γ)), ψ = ψ)

τsol = zeros(cntry_prm.Ncntry)

θ = make_θ(cntry, Rsol[cntry], Wsol[cntry], p, τsol[cntry], foo_hh_prm; points = 3, order = 1)

ω = make_ω(hh[cntry], dist[cntry], L[cntry], p, foo_hh_prm)
# makes the expenditure weights

agθ = aggregate_θ(θ, ω, cntry, foo_hh_prm)

deleteat!(agθ, cntry)

deleteat!(p, cntry)

plot(p, -agθ, seriestype = :scatter, alpha = 0.75,
    xlabel = "price",
    ylabel = "elasticity",
    legend = false)

cntrytrade = tradeshare[cntry,:]

deleteat!(cntrytrade, cntry)

df = DataFrame(θij = agθ,
               p = p,
               trade = cntrytrade,
               );

root = rootfile*"elasticity-by-partner-"*string(cntry)*".csv"

# ####################################################################################
# ####################################################################################
# # Let's construct bilateral trade elasticities

# CSV.write(root, df);

p = make_p(Wsol[1:end], TFP, d[cntry, :], cntry_prm.tariff[cntry, :] )
# prices from the perspective of those in that country

mpc = make_mpc(hh[cntry], Rsol[cntry], Wsol[cntry], p, 0.016/2, foo_hh_prm)

τeqv = zeros(foo_hh_prm.Na, foo_hh_prm.Nshocks);

fooX = make_Xsection(Rsol[cntry], Wsol[cntry], p, hh[cntry], dist[cntry],
         θ, mpc, τeqv, cntry, foo_hh_prm; Nsims = 100000);

microm = cal_make_stats(fooX, prctile = [20.0, 80.0])