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

grvdata = gravity(dftrade, display = true);

trc = trade_costs(grvdata.dist_coef, grvdata.lang_coef, grvdata.θm)

################################################################

initial_x = vec([-1.0*ones(18); zeros(18) ; -ones(6); trc.lang_coef])

gravity_as_guide(initial_x, grvdata, grv_params)

f(x) = gravity_as_guide(x, grvdata, grv_params);

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

################################################################
Ncntry = 19

T = [exp.(sol.x[1:(Ncntry - 1)]); 1] # S's are normalized -> only have 18 degrees of freedom on Ts
    
θm = [sol.x[Ncntry:((Ncntry - 1)*2)]; -sum(sol.x[Ncntry:((Ncntry - 1)*2)])] # same with this, they sum zero

dist_coef = sol.x[((Ncntry - 1)*2 + 1):((Ncntry - 1)*2 + 6)] # six distance bins

lang_coef = sol.x[((Ncntry - 1)*2 + 7):end] # the language stuff


trc = trade_costs(dist_coef, lang_coef, θm)

grv, W, πshares, dfmodel = gravity_as_guide(trc, T, grv_params.dfcntryfix, grv_params, solver = false)

################################################################


plot(dfmodel.trade, dftrade.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)

