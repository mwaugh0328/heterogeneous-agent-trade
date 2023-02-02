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

trc = trade_costs(grvdata.dist_coef, grvdata.lang_coef, grvdata.θm) 

# ################################################################
# # Recover the trade costs and technology parameters

θ = 4.0

d = zeros(19,19)
T = zeros(19)
W = ones(19)

# make_trade_costs!(dfcntryfix, trc, d, θ)

 make_technology!(grvdata, T, W, θ)

# W = solve_trade_balance(dflabor.L, d, T, θ)

# # πshares, Φ = eaton_kortum(W, d, T, θ)

initial_x = vec([-1.0*ones(18); zeros(18) ; -ones(6); trc.lang_coef])

gravity_as_guide(initial_x, dfcntryfix, dflabor.L, grvdata)

f(x) = gravity_as_guide(x, dfcntryfix, dflabor.L, grvdata);

function f!(fvec, x)

    fvec .= f(x)

end

n = length(initial_x)
diag_adjust = n - 1

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

print(sol)

# # ################################################################
# # ################################################################

# test = gravity_as_guide(trc, W, T, dfcntryfix, dflabor.L)

# W = zeros(18)
# T = log.(T)

# initial_x = vec([T; W; trc.θm; trc.dist_coef; trc.lang_coef])

# foo = gravity_as_guide(initial_x, dfcntryfix, dflabor.L, grvdata)

# f(x) = gravity_as_guide(x, dfcntryfix, dflabor.L, grvdata);

# function f!(fvec, x)

#     fvec .= f(x)

# end

# ###################################################################

# n = length(initial_x)
# diag_adjust = n - 1

# sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
#       ml=diag_adjust, mu=diag_adjust,
#       diag=ones(n),
#       mode= 1,
#       tol=1e-5,
#        )

# print(sol)
