include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("static-trade-environment.jl")
using MINPACK

Ncntry = 4

dtest = 2.0
d = dtest.*ones(Ncntry,Ncntry)
d[diagind(d)] .= 1.0

TFP = [1.0, 1.0, 0.75, 1.25]

mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, Nshocks = 5, 
γ = 2.0, ϕ = 3, amax = 8.0, σ = 0.40, ρ = 0.20, σϵ = 0.25, d = d, TFP = TFP)

@unpack Na, Nshocks, Ncntry, TFP = mdl_prm

R = 1.03*ones(Ncntry);
W = 1.0*ones(Ncntry);

f(x) = world_equillibrium(x, mdl_prm, hh_solution_method = "itteration");

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

Y, tradeflows, A_demand, tradeshare, hh, dist = world_equillibrium(Rsol, Wsol, mdl_prm, hh_solution_method = "itteration");