include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("static-trade-environment.jl")
using Plots



####################################################################################
####################################################################################
# Ncntry = 19

# dftrade = DataFrame(CSV.File("ek-trade.csv"))

# d = reshape(dftrade.d, Ncntry,Ncntry)

# df = DataFrame(CSV.File("solution-fg.csv"))

# initial_x = [df.wage[2:end]; 1.035]

TFP = [10.10; 1.0]
wage = [10.10; 1.0]

d_ij = 15.5
d = [1.0 d_ij; d_ij 1.0]

Ncntry = size(d)[1]

hh_prm = household_params(Ncntry = 2, Na = 100, 
γ = 1.5, ϕ = 0.5, amax = 8.0, σϵ = 0.25, β = 0.92)

agrid = make_agrid(hh_prm, TFP[1])

foo = household_params(hh_prm, agrid = agrid, TFP = TFP[1])

#the way this grid is setup seems to work

p = (wage[1:end] ./ TFP).*d[1,:]

hh = solve_household_problem(1.00, wage[1], p, foo)

dist = make_stationary_distribution(hh, foo)

adist = get_distribution(dist.state_index, dist.λ);

plot(foo.agrid , adist, alpha = 0.5, lw = 4,
    color = "dark blue", ylabel = "Probability Mass", 
    xlabel = "Asset Holdings / Avg. Income", label = false)



# @unpack Na, Nshocks, Ncntry, β, σϵ, TFP = foo

# gc = repeat(range(0.1,3,Na),1,Nshocks,Ncntry)
