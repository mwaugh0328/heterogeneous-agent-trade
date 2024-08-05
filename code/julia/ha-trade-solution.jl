struct household{T}
    asset_policy::Array{T} # asset_policy
    cons_policy::Array{T} # asset_policy
    πprob::Array{T} # choice probabilities
    Tv::Array{T} # value function
end


struct distribution{T}
    Q::Array{T} # transition matrix
    λ::Array{T} # λ
    state_index::Array{Tuple{Int64, Int64}} # index of states, lines up with λ
end


# ##########################################################################
# functions used to find a solution
# 

function world_equillibrium(x, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration")

    @unpack Ncntry = cntry_params

    @assert length(x) ≈ ( Ncntry + (Ncntry - one(Ncntry)) )

    W = [x[1:(Ncntry - one(Ncntry))]; 1.0 ]
    W = W ./ ( sum(W) / Ncntry)
    
    R = x[Ncntry:end]

    Y, tradeflows, A_demand = world_equillibrium(R, W, hh_params, cntry_params; tol_vfi = tol_vfi, tol_dis = tol_dis, 
        hh_solution_method = hh_solution_method, stdist_sol_method=stdist_sol_method)[1:3]

    goods_market = Y .- vec(sum(tradeflows, dims = 1))
    # so output (in value terms) minus stuff being purchased by others (value terms so trade costs)
    # per line ~ 70 below, if we sum down a row this is the world demand of a countries commodity. 

    asset_market = A_demand

    return [asset_market; goods_market[2:end]]

end

# ##########################################################################
# ##########################################################################
function world_equillibrium_FG(x, R, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration")
    # multiple dispatch version to find β, given R, so that asset market clears

    β = x[end]

    foo_hh_params = household_params(hh_params, β = β)

    x[end] = R # this is bc x is setup for R to be there 
    # should change later

    return world_equillibrium_FG(x, foo_hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
        hh_solution_method = "itteration", stdist_sol_method = "itteration")

end



function world_equillibrium_FG(x, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration")

    @unpack Ncntry = cntry_params

    @assert length(x) ≈ Ncntry
    

    W = [x[1:(Ncntry - one(Ncntry))]; 1.0 ]
    W = W ./ ( sum(W) / Ncntry)
    # There is an assumption that the final country is the one we 
    # are normalizing everything to. 

    R = ones(Ncntry)*x[end]

    # dfguess = DataFrame(W = W, R = R);

    # CSV.write("current-price.csv", dfguess)

    Y, tradeflows, A_demand = world_equillibrium(R, W, hh_params, cntry_params; tol_vfi = tol_vfi, tol_dis = tol_dis, 
        hh_solution_method = hh_solution_method, stdist_sol_method=stdist_sol_method)[1:3]

    goods_market = Y .- vec(sum(tradeflows, dims = 1))
    # so output (in value terms) minus stuff being purchased by others (value terms so trade costs)
    # per line ~ 70 below, if we sum down a row this is the world demand of a countries commodity. 

    asset_market = A_demand

    return [sum(asset_market); goods_market[2:end]]

end

# ##########################################################################
# ##########################################################################

function world_equillibrium_FG_τ(x, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration")
    # for the situation in which there are rebates to households

    @unpack Ncntry = cntry_params

    #@assert length(x) ≈ Ncntry

    W, τ, R = unpack_xvec(x, Ncntry)

    #dfguess = DataFrame(W = W, R = R);

    #CSV.write("current-price.csv", dfguess)

    Y, tradeflows, A_demand, Gbudget = world_equillibrium(R, W, τ, hh_params, cntry_params; tol_vfi = tol_vfi, tol_dis = tol_dis, 
        hh_solution_method = hh_solution_method, stdist_sol_method=stdist_sol_method)[1:4]

    goods_market = Y .- vec(sum(tradeflows, dims = 1))
    # so output (in value terms) minus stuff being purchased by others (value terms so trade costs)
    # per line ~ 70 below, if we sum down a row this is the world demand of a countries commodity. 

    asset_market = A_demand

    return [sum(asset_market); goods_market[2:end]; Gbudget]

end


# ##########################################################################
# ##########################################################################

function world_equillibrium(R, W, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration")
    #multiple dispatch here for the case with no transfer

    τ = zeros(cntry_params.Ncntry)

    return world_equillibrium(R, W, τ, hh_params, cntry_params, tol_vfi = tol_vfi, tol_dis = tol_dis, 
        hh_solution_method = hh_solution_method, stdist_sol_method = stdist_sol_method)    

end

# ##########################################################################
# ##########################################################################

function world_equillibrium(R, W, τ, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration")

    @assert hh_params.ϕ > 0.0

    @unpack Ncntry, TFP, d, tariff, L = cntry_params
    @unpack ψslope, γ, σϵ = hh_params

    @assert length(R) ≈ Ncntry
    @assert length(W) ≈ Ncntry

    Y = similar(W)
    A_demand = similar(R)
    Gbudget = similar(R)

    tradeflows = Array{Float64}(undef,Ncntry,Ncntry)
    tradeshare = Array{Float64}(undef,Ncntry,Ncntry)

    hh = Array{household{Float64}}(undef,Ncntry)
    dist = Array{distribution{Float64}}(undef,Ncntry)

    Threads.@threads for cntry = 1:Ncntry

        p = make_p(W, TFP, d[cntry, :], tariff[cntry, :] )

        ψ = make_ψ(cntry, ψslope.*TFP[cntry].^(1.0 - γ), hh_params)
        # this creates the z quality shifter
        # scaled in a way that is invariant to level of TFP

        agrid = make_agrid(hh_params, TFP[cntry])
        # this creates teh asset grid so it's alwasy a fraction of home labor income

        foo_hh_params = household_params(hh_params, agrid = agrid, 
                TFP = TFP[cntry], L = L[cntry], σϵ = σϵ*(TFP[cntry]^(1.0 - γ)), ψ = ψ)

        hh[cntry], dist[cntry] = compute_eq(R[cntry], W[cntry], p, τ[cntry], foo_hh_params, tol_vfi = tol_vfi, tol_dis = tol_dis,
            hh_solution_method = hh_solution_method, stdist_sol_method = stdist_sol_method)

    end

    for cntry = 1:Ncntry

        p = make_p(W, TFP, d[cntry, :], tariff[cntry, :] )

        ψ = make_ψ(cntry, ψslope.*TFP[cntry].^(1.0 - γ), hh_params)

        agrid = make_agrid(hh_params, TFP[cntry])

        foo_hh_params = household_params(hh_params, agrid = agrid, 
                TFP = TFP[cntry], L = L[cntry], σϵ = σϵ*(TFP[cntry]^(1.0 - γ)), ψ = ψ)

        output, tradestats = aggregate(R[cntry], W[cntry], p, τ[cntry], tariff, cntry, hh[cntry], dist[cntry], foo_hh_params)

        Y[cntry] = output.production

        tradeflows[cntry, :] = tradestats.bilateral_imports

        tradeshare[cntry, :] = tradestats.bilateral_imports ./ output.PC
        # the way I read this is fix a row, then across the columns this is how much cntry in position cntry
        # is buying/importing from of the other commodities. 

        Gbudget[cntry] = tradestats.tariff_revenue - output.G 

        A_demand[cntry] = output.Aprime

    end

return Y, tradeflows, A_demand, Gbudget, tradeshare, hh, dist

end

# ##########################################################################
# ##########################################################################
function make_agrid(hh_params, TFP)

    return convert(Array{Float64, 1}, range(-hh_params.ϕ*TFP, hh_params.amax*TFP, length = hh_params.Na))

end


function micro_trade_elasticity(R, W, p, τ, home, source, model_params; tol_vfi = 1e-6, hh_solution_method = "itteration")

    hh = solve_household_problem(R, W, p, τ, model_params, tol = tol_vfi, solution_method = hh_solution_method)

    return (hh.πprob[:,:,source] ./ hh.πprob[:,:,home])

end

function micro_trade_elasticity(W, p, home, source, model_params)
    #using mulitiple dispatch here... if no R, then it's 
    # the static model

    πprob = logit_trade(W, p, model_params)[2]

    return (πprob[:,:,source] ./ πprob[:,:,home])

end

# ##########################################################################
# ##########################################################################

function compute_eq(R, W, p, τ, model_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration")
    # Does everything...
    # (1) Sovles hh problem
    # (2) Constructs stationary distribution

    hh = solve_household_problem(R, W, p, τ, model_params, tol = tol_vfi, solution_method = hh_solution_method)

    #println("hh problem solved")

    dist = make_stationary_distribution(hh, model_params, tol = tol_dis, solution_method = stdist_sol_method)

    #println("stationatry distribution found")
    
return hh, dist

end

# ##########################################################################
# ##########################################################################

function make_stationary_distribution(household, model_params; tol = 1e-10, solution_method = "itteration") 

@unpack Na, Nshocks = model_params

Q = Array{Float64}(undef, Na*Nshocks, Na*Nshocks)

make_Q!(Q, household, model_params)

if solution_method == "nl-fixedpoint"

    g(λ) = law_of_motion(λ, transpose(Q))

    initialvalue = zeros(size(Q)[1], 1)

    initialvalue .= 1.0 / size(Q)[1]

    solution = fixedpoint(g, initialvalue, ftol = tol, method = :anderson)

    λ = solution.zero

elseif solution_method == "itteration"

    λ = itterate_stationary_distribution(Q; tol = tol)

end

state_index = Array{Tuple{eltype(Na), eltype(Na)}}(undef, Na*Nshocks, 1)

make_state_index!(state_index, model_params)

return distribution(Q, λ, state_index)

end

##########################################################################
##########################################################################
function itterate_stationary_distribution(Q; tol = 1e-10, Niter = 5000) 
    # this is faster than the quant econ canned routine
    # from lyon-waugh implementation. Not as fast as using
    # NLsolve fixedpoint routine below.
    
    # Takes the transition matrix above then we know a stationary distribution
    # must satisfy the fixed point relationsip λ = Q'*λ  
    
    λ = zeros(size(Q)[1], 1)
    λ = convert(Array{eltype(Q)}, λ)
    
    λ .= 1.0 / size(Q)[1]
    
    Lnew = similar(λ)
    Qtran = transpose(Q)
    
    for iter in 1:Niter
        
        #Lnew = transpose(Q) * λ
        
        #Lnew = law_of_motion(λ, transpose(Q))

        mul!(Lnew, Qtran, λ) 
        #this ordering is also better in julia
        # than my matlab implementation of Q*λ (1, na*nshock)
                
        err = maximum(abs, λ - Lnew)
        
        copy!(λ, Lnew)
        # this surprisingly makes a big difference
        # but in the vfi it causes a slowdown?
        
        err < tol && break
        
        if iter == Niter

            println("distribution may not have converged")
            println("check the situation")
    
        end
        
    end

    return λ

end

# ##########################################################################

function solve_household_problem(R, W, p, τ, model_params; tol = 10^-6, solution_method = "itteration")

    if solution_method == "nl-fixedpoint"

        Kga, Kgc, πprob, Tv = policy_function_fixedpoint(R, W, p, τ, model_params; tol = tol)

    elseif solution_method == "itteration"

        Kga, Kgc, πprob, Tv = policy_function_itteration(R, W, p, τ, model_params; tol = tol, Niter = 500)

    elseif solution_method == "itteration_old"

        Kga, Kgc, πprob, Tv = policy_function_itteration_old(R, W, p, τ, model_params; tol = tol, Niter = 500)
        
    end
    
    return household(Kga, Kgc, πprob, Tv)

end

##########################################################################
##########################################################################

# This function is renamed to '_old' just to check with old codes,
# will delete once all test is done
function policy_function_itteration_old(R, W, p, τ, model_params; tol = 10^-6, Niter = 500)
    
    @unpack Na, Nshocks, Ncntry, β, σϵ, ψ = model_params

    # this is the guess... always start at borrowing cosntraint
    gc = Array{Float64}(undef, Na, Nshocks, Ncntry)
    Kgc = similar(gc)
    
    make_gc_guess!(gc, R, W, p, model_params)
    
    #gc = repeat(range(0.1,3,Na),1,Nshocks,Ncntry)
    #println(gc)

    v = -ones(Na, Nshocks, Ncntry)/(1-β)
    Tv = similar(v)

    for iter in 1:Niter
        
        Kgc, Tv = coleman_operator_old(gc, v, R, W, p, τ, model_params)[1:2]

        err = vec_max(Kgc, gc)

        #println(iter)

        if err < tol

            #println(iter)

            break
        end

        copy!(gc, Kgc)

        copy!(v,Tv)

        if iter == Niter

          println("value function may not have converged")
          println("check the situation")
          
        end

    end

    Kgc, Tv, Kga = coleman_operator_old(gc, Tv, R, W, p, τ, model_params)

    πprob = make_πprob(Tv, σϵ, ψ)

    return Kga, Kgc, πprob, Tv
    
end

# We don't change this function at all as we have 'coleman_operator' multiple dispatch
function policy_function_itteration(R, W, p, τ, model_params; tol = 10^-6, Niter = 500)
    
    @unpack Na, Nshocks, Ncntry, β, σϵ, ψ = model_params

    # this is the guess... always start at borrowing cosntraint
    gc = Array{Float64}(undef, Na, Nshocks, Ncntry)
    Kgc = similar(gc)
    
    make_gc_guess!(gc, R, W, p, model_params)
    
    #gc = repeat(range(0.1,3,Na),1,Nshocks,Ncntry)
    #println(gc)

    v = -ones(Na, Nshocks, Ncntry)/(1-β)
    Tv = similar(v)

    for iter in 1:Niter
        
        Kgc, Tv = coleman_operator(gc, v, R, W, p, τ, model_params)[1:2]

        err = vec_max(Kgc, gc)

        #println(iter)

        if err < tol

            #println(iter)

            break
        end

        copy!(gc, Kgc)

        copy!(v,Tv)

        if iter == Niter

          println("value function may not have converged")
          println("check the situation")
          
        end

    end

    Kgc, Tv, Kga = coleman_operator(gc, Tv, R, W, p, τ, model_params)

    πprob = make_πprob(Tv, σϵ, ψ)

    return Kga, Kgc, πprob, Tv
    
end

function vec_max(x,y)
    
    s = 0.0

    @inbounds @turbo for i ∈ eachindex(x)

        s = ifelse(abs(x[i]-y[i]) > s, abs(x[i]-y[i]), s)

    end
    s
end

##########################################################################
##########################################################################

# this function needs double check
function policy_function_fixedpoint(R, W, p, τ, model_params; tol = 10^-6)

    @unpack Na, Nshocks, Ncntry, statesize, β, σϵ, ψ, TFP = model_params

    #foo, foobar = policy_function_itteration(R, W, p, model_params, Niter = 2)
    #policy_o = vcat(foo, foobar)

    Kgc = repeat(range(0.1,3,Na)*TFP,1,Nshocks,Ncntry)

    println(Kgc)

    Tv = -ones(Na, Nshocks, Ncntry)/(1-β)

    policy_o = vcat(Kgc, Tv)

    # have seen some convergence issues sometimes

    K(policy) = coleman_operator(policy, R, W, p, τ, model_params)

    solution = fixedpoint(K, policy_o, ftol = tol, method = :anderson);
    
    if solution.f_converged == false
        println("did not converge")
    end

    #Kgc = solution.zero[1:Na, :, :]

    #Tv = solution.zero[(Na+1):end, :, :]

    Kgc, Tv, Kga = coleman_operator(solution.zero[1:Na, :, :], solution.zero[(Na+1):end, :, :], R, W, p, τ,  model_params)[2:3]

    πprob = make_πprob(Tv, σϵ, ψ)

    return Kga, Kgc, πprob, Tv

end

##########################################################################
##########################################################################


function make_gc_guess!(gc, R, W, p, model_params)

    @unpack mc, Ncntry, Nshocks, agrid = model_params

    shocks = exp.(mc.state_values)
    
    for cntry = 1:Ncntry

        for shk = 1:Nshocks

            gc[:, shk, cntry] .= (-agrid[1]  + (R*agrid[1] + W*shocks[shk]) ) / p[cntry]
       
        end

    end

end

##########################################################################
##########################################################################

function value_function_fixedpoint(R, W, p, τ, model_params; tol = 10^-6)

    @unpack Ncntry, Na, Nshocks, β, mc, σϵ = model_params

    u = Array{eltype(R)}(undef, Na, Na, Nshocks, Ncntry)
    
    make_utility!(u, R, W, p, τ, model_params)

# define the inline function on the bellman operator. 
# so the input is v (other stuff is fixed)
    
    T(v) = bellman_operator_upwind(v, u, mc.p, β, σϵ)

    Vo = -ones(Na, Nshocks) / (1-β); #initial value
    Vo = convert(Array{eltype(u)}, Vo)

    solution = fixedpoint(T, Vo, ftol = tol, method = :anderson);
    
    if solution.f_converged == false
        println("did not converge")
    end

    πprob = make_πprob(solution.zero, σϵ)

    return solution.zero, πprob, u

end







