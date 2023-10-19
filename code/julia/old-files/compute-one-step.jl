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
#################################################################################

df = DataFrame(CSV.File("solution-fg.csv"))

L = df.L

Ncntry = size(L)[1]

γ = 1.5

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = γ, ϕ = 0.5, amax = 8.0, σϵ = 0.25)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

dfparams = DataFrame(CSV.File("current-guess-ek-new.csv"))
#dfparams = DataFrame(CSV.File("current-guess-log-ek.csv"))

initial_x = dfparams.guess

# initial_x = log.([TFP[1:18]; 1.02])

foo  = [exp.(initial_x[1:18]) ; 1.02; initial_x]

calibrate_world_equillibrium(foo, grvdata, grv_params, hh_prm, cntry_prm )


f(x) = calibrate_world_equillibrium(x, grvdata, grv_params, hh_prm, cntry_prm )

function f!(fvec, x)

    fvec .= f(x)

end

n = length(foo)

diag_adjust = n - 1

sol = fsolve(f!, foo, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-3,
       )
