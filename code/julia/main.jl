include("ha-trade.jl")

γ = 1.5 # curvatuve on CRRA utility function
σϵ = 0.25 # logit dispersion parameter
Ncntry = 2 # number of countries

# this setups up parameters on the household side
hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92, γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ);

Na = hh_prm.Na
Nshocks = hh_prm.Nshocks
β = hh_prm.β

TFP = [1.0; 1.0]

τ = [0.0; 0.0]

L = [1.0; 1.0]

d_ij = 1.745

d = [1.0 d_ij; d_ij 1.0]

init_c = Array{Float64}(undef, Na, Nshocks, Ncntry)
init_R = 1.0
init_W = 1.02
init_p = make_p(init_W, TFP, d[1, :], τ[1, :])
make_gc_guess!(init_c, init_R, init_W, init_p, hh_prm)
init_v = -ones(Na, Nshocks, Ncntry)/(1-β)


coleman_operator(init_c, init_v, init_R, init_W, init_p, τ[1], hh_prm)