include("ha-trade.jl")
using Plots

####################################################################################
####################################################################################

TFP = [10.0; 1.0]
wage = [10.0; 1.0]

d_ij = 10.5
d = [1.0 d_ij; d_ij 1.0]

Ncntry = size(d)[1]

γ = 1.5
σϵ = 0.25

hh_prm = household_params(Ncntry = 2, Na = 100, 
γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ, β = 0.92, ρ = 0.90)

agrid = make_agrid(hh_prm, TFP[1])

foo = household_params(hh_prm, agrid = agrid, TFP = TFP[1], σϵ = σϵ*(TFP[1]^(1.0 - γ)))

#the way this grid is setup seems to work

p = (wage[1:end] ./ TFP).*d[1,:]
R = 1.00

τ = 0.0

@time hh = solve_household_problem(R, wage[1], p, τ, foo)

@time dist = make_stationary_distribution(hh, foo)

θelas = make_θ(1, R, wage[1], p, τ, foo)

adist = get_distribution(dist.state_index, dist.λ);

mpc = make_mpc(hh, R, wage[1], p, τ + 0.1, foo)

plot(foo.agrid , adist, alpha = 0.5, lw = 4,
    color = "dark blue", ylabel = "Probability Mass", 
    xlabel = "Asset Holdings / Avg. Income", label = false)




