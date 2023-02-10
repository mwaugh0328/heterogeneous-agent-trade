include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("static-trade-environment.jl")

using FiniteDifferences

using Plots

Ncntry = 5

dtest = 1.25
d = dtest.*ones(Ncntry,Ncntry)
d[diagind(d)] .= 1.0

mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, Nshocks = 5, 
γ = 2.0, ϕ = 3, amax = 8.0, σ = 0.10, ρ = 0.95, σϵ = 0.25, d = d)

@unpack Na, Nshocks, Ncntry, TFP = mdl_prm


R = 1.03*ones(Ncntry);
W = 1.0*ones(Ncntry);

cntry = 1

p = (W[cntry] ./ TFP[cntry] ) .* d[cntry, :]

gc, πprob, v = logit_trade(W[cntry], p, mdl_prm);

# g(x) = make_Mij_Mii(x[1], W[cntry], cntry, 2, mdl_prm)

# ∂M_∂d = jacobian(central_fdm(2, 1), g, [1.5])[1]

# ∂logM_∂logd = reshape(∂M_∂d, Na, Nshocks) .* (1.5 ./ g(1.5))




function make_Mij_Mii(∂d, W, home, source, model_params)

    p = make_p(∂d, W, home, source, model_params)

    gc, πprob = logit_trade(W[home], p, model_params)[1:2]

    Mij_Mii = (πprob[:,:,source].*gc[:,:,source].*p[source]) ./ (πprob[:,:,home].*gc[:,:,home].*p[home])

    return Mij_Mii

end

function make_Mij_Mii(∂d, R, W, home, source, model_params)

    p = make_p(∂d, W, home, source, model_params)

    hh = solve_household_problem(R, W, p, model_params, solution_method = "nl-fixedpoint", tol = 1e-10)
    #itteration is not working properly...need to fix, nlfixed point yes

    state_index = Array{Tuple{eltype(Int64), eltype(Int64)}}(undef, model_params.Na*model_params.Nshocks, 1)

    make_state_index!(state_index, model_params)

    gc, πprob = get_trade(R, W, hh.asset_policy, hh.πprob, state_index, model_params)

    #Mij_Mii = (πprob[:,home].*gc[:,source]) ./ (πprob[:,source].*gc[:,home])

    Mij_Mii = (gc[:,source]) ./ (gc[:,home])

    return Mij_Mii

end


# function make_Mij_Mii(∂d, R, W, home, source, model_params)

#     p = make_p(∂d, W, home, source, model_params)

#     hh = solve_household_problem(R[home], W[home], p, mdl_prm)

#     Mij_Mii = (gc[:,:,source].*p[source]) ./ (gc[:,:,home].*p[home])

#     return Mij_Mii

# end

function make_p(∂d, W, home, source, model_params)

    @unpack TFP, d = model_params

    p = (W ./ TFP[home] ) .* d[home, :]

    p[source] = (W ./ TFP[home] )*∂d

    return p

end


h(x) = make_Mij_Mii(x[1], R[cntry], W[cntry], cntry, 2, mdl_prm)

∂M_∂d = jacobian(central_fdm(4, 1), h, [dtest])[1]

∂logM_∂logd = ∂M_∂d .* (dtest ./ h(dtest))



# hh = solve_household_problem(R[cntry], W[cntry], p, mdl_prm)

# state_index = Array{Tuple{eltype(Na), eltype(Na)}}(undef, Na*Nshocks, 1)

# make_state_index!(state_index, mdl_prm)

# c, share = get_trade(R[cntry], W[cntry], hh.asset_policy, hh.πprob, state_index, mdl_prm)


# f(x) = Mij_Mii(W, x[1], 2, 1, mdl_prm)

# ∂M_∂d = jacobian(central_fdm(2, 1), f, [d])[1]

# trade_elasticity = reshape(∂M_∂d, Na, Nshocks) .* (d ./ f(d))
