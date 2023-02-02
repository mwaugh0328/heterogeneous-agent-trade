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

# ################################################################
# reconstruct the output 

T = [exp.(sol.x[1:18]); 1] 
θm = [sol.x[19:(36)]; -sum(sol.x[19:(36)])] 
dist_coef = sol.x[37:37+5]
lang_coef = sol.x[43:end]

trc = trade_costs(dist_coef, lang_coef, θm)

grv, W, πshares = gravity_as_guide(trc, T, dfcntryfix, dflabor.L, 4.0, solver = false)

# # ################################################################


