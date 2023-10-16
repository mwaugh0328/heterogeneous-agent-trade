include("ha-trade.jl")
using MINPACK
using Plots
using CSV
using DataFrames
using StatsBase

####################################################################################
####################################################################################
# This sets up the EK trade data and gravity stuff

trade_cost_type = "ek"

grvdata, grv_params, L, dftrade = make_gravity_params(trade_cost_type)

####################################################################################
#################################################################################

Ncntry = size(L)[1]

γ = 1.45
σϵ = 0.33
ψslope = 0.725

# γ = 1.00
# σϵ = 0.2369
# ψslope = 0.0

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ, ψslope = ψslope)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

R = 1.01

dfparams = DataFrame(CSV.File("./calibration-files/current-guess-145-30-70.csv"))
# dfparams = DataFrame(CSV.File("./calibration-files/current-guess-log-22.csv"))

dfWguess = DataFrame(CSV.File("./calibration-files/current-wage-guess-145-30-70.csv"))
# dfWguess = DataFrame(CSV.File("./calibration-files/current-wage-guess-log-22.csv"))

micro_moment = 0.0

initial_x = [dfWguess.guess[1:19]; dfparams.guess[1:end]]
# last one in the W guess is the discount factor

out = calibrate_world_equillibrium(initial_x, R, grvdata, grv_params, hh_prm, 
                                cntry_prm, trade_cost_type = trade_cost_type)


f(x) = calibrate_world_equillibrium(x, R, grvdata, grv_params, hh_prm, 
                cntry_prm, trade_cost_type = trade_cost_type)

function f!(fvec, x)

    fvec .= f(x)

end

n = length(initial_x)

diag_adjust = 19

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,)



# dfguess = DataFrame(guess = sol.x[20:end]);
# CSV.write("./calibration-files/current-guess-log-24.csv", dfguess)

# dfWRguess = DataFrame(guess = sol.x[1:19]);
# CSV.write("./calibration-files/current-wage-guess-log-24.csv", dfWRguess)



