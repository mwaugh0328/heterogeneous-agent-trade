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


trade_cost_type = "waugh"

grvdata = gravity(dftrade, display = true, trade_cost_type = trade_cost_type );

trc = trade_costs(grvdata.dist_coef, grvdata.lang_coef, grvdata.θm)

grv_params = gravity_params(L = dflabor.L, dfcntryfix = dfcntryfix, Ncntry = 19)

####################################################################################
#################################################################################

df = DataFrame(CSV.File("solution-fg.csv"))

L = df.L

Ncntry = size(L)[1]

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = 1.5, ϕ = 0.5, amax = 8.0, σϵ = 0.25)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

initial_x = [vec(log.(df.TFP[1:18])); trc.θm[1:18] ; trc.dist_coef; trc.lang_coef]

out, Wsol, Rsol, πshare = calibrate(initial_x, grvdata, grv_params, hh_prm, 
                                cntry_prm, trade_cost_type = trade_cost_type) 

price_guess = log.([Wsol[1:18]; Rsol[1]])

f(x) = calibrate(x, grvdata, grv_params, hh_prm, cntry_prm, trade_cost_type = trade_cost_type)[1]

function f!(fvec, x)

    fvec .= f(x)

end

n = length(initial_x)

diag_adjust = n - 1

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-3,
       )



