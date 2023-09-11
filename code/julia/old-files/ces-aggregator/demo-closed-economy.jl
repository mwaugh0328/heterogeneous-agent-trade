include("ha-trade.jl")

Pces = 1.0
W = 1.0
τ_rev = 0.0
R = 1.029

nolabor = model_params(ϑ = 0.0, σw = 0.001, σa = 0.005)

hh, dist = compute_eq(Pces, W, τ_rev, R, nolabor)