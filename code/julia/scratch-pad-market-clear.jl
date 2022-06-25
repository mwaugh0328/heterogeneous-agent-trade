include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("static-trade-environment.jl")
using MINPACK
using Plots
using CSV
using DataFrames

####################################################################################
# brings in EK data 
Ncntry = 19

dftrade = DataFrame(CSV.File("ek-trade.csv"))

d = reshape(dftrade.d, Ncntry,Ncntry)

df = DataFrame(CSV.File("solution.csv"))

initial_x = [df.wage[2:end]; df.interest_rate]

#initial_x =[df.TFP[2:end]; 1.02*ones(Ncntry)]

TFP = df.TFP
L = df.L

####################################################################################

Ncntry = size(d)[1]


mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, 
γ = 1.5, ϕ = 2.0, amax = 8.0, σϵ = 0.25, d = d, TFP = TFP, L = L)

@unpack Na, Nshocks, Ncntry, TFP = mdl_prm

# R = 1.022*ones(Ncntry);
# W = TFP;

f(x) = world_equillibrium(x, mdl_prm, hh_solution_method = "itteration", stdist_sol_method = "itteration");

function f!(fvec, x)

    fvec .= f(x)

end

###################################################################

n = length(initial_x)
diag_adjust = n - 1

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

print(sol)

###################################################################

Wsol = [1.0; sol.x[1:Ncntry-1]]
Rsol = sol.x[Ncntry:end]

 Y, tradeflows, A_demand, tradeshare, hh, dist = world_equillibrium(Rsol, Wsol, 
     mdl_prm, hh_solution_method = "itteration");

plot(log.(vec(tradeshare)), log.(dftrade.tradesharedata), seriestype = :scatter)


df = DataFrame(country_index = range(1,19), 
     wage = Wsol,
     interest_rate = Rsol,
     TFP = TFP,
     L = L
     );

CSV.write("solution.csv", df)

adist = get_distribution(dist[19].state_index, dist[19].λ);

plot(mdl_prm.agrid, adist, alpha = 0.5, lw = 4,
    color = "dark blue", ylabel = "Probability Mass", xlabel = "Asset Holdings", label = false)