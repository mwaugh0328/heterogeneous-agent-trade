include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("static-trade-environment.jl")
using MINPACK

Ncntry = 4

dtest = 1.75
d = dtest.*ones(Ncntry,Ncntry)
d[diagind(d)] .= 1.0

mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, Nshocks = 5, 
γ = 1.5, ϕ = 3, amax = 8.0, σ = 0.40, ρ = 0.20, σϵ = 0.25, d = d)

@unpack Na, Nshocks, Ncntry, TFP = mdl_prm

R = 1.037*ones(Ncntry);
W = 0.50*ones(Ncntry);

f(x) = world_equillibrium(x, mdl_prm, hh_solution_method = "itteration");

function f!(fvec, x)

    fvec .= f(x)

end

initial_x = [W; R]

n = length(initial_x)
diag_adjust = n - 1

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )