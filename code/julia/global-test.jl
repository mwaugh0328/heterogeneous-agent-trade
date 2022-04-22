include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")

Ncntry = 12
d = 1.5.*ones(Ncntry,Ncntry)
d[diagind(d)] .= 1.0

mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, Nshocks = 5, 
γ = 3.0, ϕ = 3, amax = 8.0, σ = 0.3919, ρ = 0.20, d = d)

@unpack Ncntry = mdl_prm

R = 1.029*ones(Ncntry);
W = 1.0*ones(Ncntry);

@time Y, tradeflows, Ademand = world_equillibrium(R, W, mdl_prm)


# @time Y, tradeflows, Ademand = world_equillibrium(R, W, mdl_prm)
#   0.928947 seconds (9.30 M allocations: 1.170 GiB, 14.82% gc time)