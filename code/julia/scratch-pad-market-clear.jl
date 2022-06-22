include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("static-trade-environment.jl")
using MINPACK
using MAT
using Plots

####################################################################################
# brings in EK data 
file = matopen("../../ek-data/ek-output.mat")
tradesharedata = read(file, "tssdmat")'
# I set my model up so row is an importer, column is exporter and 
# that accross columns should sum to one this is the oppisite so flip

d = read(file, "rtausd")'
# same deal

L = [0.054, 0.024, 0.029, 0.094, 0.017, 0.019,
    0.181, 0.0225, 0.025, 0.159, 0.544, 0.043, 0.010, 0.015,
    0.026, 0.10, 0.031, 0.186, 1.0]

TFP = [0.36,0.30,0.22,0.47,0.32,0.41,0.61,0.75,0.14,
    0.57,0.97,0.28,0.22,0.37,0.13,0.33,0.47,0.53,1.0]

TFP .=  TFP.^(1. / 3.6)

####################################################################################

Ncntry = size(d)[1]

#TFP = vec(exp.(read(file, "ssd")))


# dtest = 2.0
# d = dtest.*ones(Ncntry,Ncntry)
# d[diagind(d)] .= 1.0

#TFP = ones(Ncntry)

mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, Nshocks = 5, 
γ = 2.0, ϕ = 3, amax = 8.0, σ = 0.15, ρ = 0.90, σϵ = 0.25, d = d, TFP = TFP, L = L)

@unpack Na, Nshocks, Ncntry, TFP = mdl_prm

R = 1.03*ones(Ncntry);
W = TFP;

f(x) = world_equillibrium(x, mdl_prm, hh_solution_method = "itteration", stdist_sol_method = "itteration");

function f!(fvec, x)

    fvec .= f(x)

end

initial_x = [W[2:end]; R]

n = length(initial_x)
diag_adjust = n - 1

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

print(sol)

Wsol = [1.0; sol.x[1:Ncntry-1]]
Rsol = sol.x[Ncntry:end]

 Y, tradeflows, A_demand, tradeshare, hh, dist = world_equillibrium(Rsol, Wsol, 
     mdl_prm, hh_solution_method = "itteration");

     plot(log.(vec(tradeshare)), log.(vec(tradesharedata)), seriestype = :scatter)