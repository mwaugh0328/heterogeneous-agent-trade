include("ha-trade.jl")
using MINPACK
using Plots
using CSV
using DataFrames

####################################################################################
####################################################################################
# This sets up the EK trade data and gravity stuff

trade_cost_type = "ek"

grvdata, grv_params, L, dftrade = make_gravity_params(trade_cost_type)

####################################################################################
#################################################################################

Ncntry = size(L)[1]

γ = 1.450
σϵ = 0.33
ψslope = 0.725

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ, ψslope = ψslope)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

R = 1.01

dfparams = DataFrame(CSV.File("./calibration-files/current-guess-145-30-725.csv"))

initial_x = dfparams.guess

out, Wsol, β, πshare = calibrate(initial_x, R, grvdata, grv_params, hh_prm, 
                                cntry_prm, trade_cost_type = trade_cost_type)[1:4]


f(x) = calibrate(x, R, grvdata, grv_params, hh_prm, cntry_prm, trade_cost_type = trade_cost_type)[1]

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



