include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")

using Plots

mdl_prm = model_params(Ncntry = 5, Na = 100, Nshocks = 5, γ = 3.0, ϕ = 3, amax = 8.0, σ = 0.3919, ρ = 0.20);

@unpack Na, Nshocks, Ncntry, β = mdl_prm

gc = ones(Na, Nshocks, Ncntry)

v = -ones(size(gc)) / (1- β)

country = 1

p = 1.5.*ones(Ncntry)

p[country] = 1.0

R = 1.029;
W = 1.0;

@time hh, dist = compute_eq(R, W, p, mdl_prm)



@time agstats, tradestats = aggregate(R, W, p, country, hh, dist, mdl_prm, display = true)

c_by_variety, variety_share = get_trade(R, W, hh.asset_policy, hh.πprob, dist.state_index, mdl_prm)

# asset_dist = get_distribution(dist.state_index, dist.λ);

# plot(mdl_prm.agrid, asset_dist)

# get_aprime(hh.asset_policy, hh.πprob, dist.state_index, mdl_prm)

# get_consumption(R, W, hh.asset_policy, hh.πprob, dist.state_index, mdl_prm)

# country = 1

# imports, import_share, home_consumption = get_trade(R, W, country, hh.asset_policy, hh.πprob, dist.state_index, mdl_prm)


# stats = aggregate(R, W, p, country, hh, dist, model_params)
#@time foo = coleman_operator(gc, v, R, W, p, mdl_prm);

# @time Kgc, πprob, Tv = policy_function_itteration(R, W, p, mdl_prm);

# @time hh = solve_household_problem(R, W, p, mdl_prm);

# Q = Array{Float64}(undef,Na*Nshocks, Na*Nshocks)

# @time make_Q!(Q, hh, mdl_prm)



#foo = log_sum_column(Tv, prams.σϵ)

#@time altTv, u = value_function_fixedpoint(R, W, p, prams);


