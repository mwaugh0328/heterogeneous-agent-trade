struct household{T}
    asset_policy::Array{T} # asset_policy
    πprob::Array{T} # choice probabilities
    Tv::Array{T} # value function
end

struct distribution{T}
    Q::Array{T} # transition matrix
    λ::Array{T} # λ
    state_index::Array{Tuple{Int64, Int64}} # index of states, lines up with λ
end

# # ##########################################################################
# # # functions used to find a solution

function world_equillibrium(x, model_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "nl-fixedpoint", stdist_sol_method = "nl-fixedpoint")

    W = x[1:Ncntry]
    R = x[Ncntry+1:end]

    Y, tradeflows, A_demand = world_equillibrium(R, W, model_params; tol_vfi = tol_vfi, tol_dis = tol_dis, 
        hh_solution_method = hh_solution_method, stdist_sol_method=stdist_sol_method)[1:3]

    goods_market = Y .- vec(sum(tradeflows, dims = 1))
    # so output (in value terms) minus stuff being purchased by others (value terms so trade costs)
    # per line ~ 70 below, if we sum down a row this is the world demand of a countries commodity. 

    asset_market = A_demand

    return [asset_market; goods_market]

end

# ##########################################################################
# ##########################################################################


function world_equillibrium(R, W, model_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "nl-fixedpoint", stdist_sol_method = "nl-fixedpoint")

@assert model_params.ϕ > 0.0

@unpack Ncntry, TFP, d = model_params

    Y = similar(W)
    A_demand = similar(R)
    tradeflows = Array{Float64}(undef,Ncntry,Ncntry)

    hh = Array{household{Float64}}(undef,Ncntry)
    dist = Array{distribution{Float64}}(undef,Ncntry)

    Threads.@threads for cntry = 1:Ncntry

        p = (W ./ TFP) .* d[cntry, :]

        hh[cntry], dist[cntry] = compute_eq(R[cntry], W[cntry], p, model_params, tol_vfi = tol_vfi, tol_dis = tol_dis,
            hh_solution_method = hh_solution_method, stdist_sol_method = stdist_sol_method)

    end

    for cntry = 1:Ncntry

        p = (W ./ TFP) .* d[cntry, :]

        output, tradestats = aggregate(R[cntry], W[cntry], p, cntry, hh[cntry], dist[cntry], model_params)

        Y[cntry] = output.production

        tradeflows[cntry, :] = tradestats.bilateral_imports
        # the way I read this is fix a row, then across the columns this is how much cntry in position cntry
        # is buying/importing from of the other commodities. 

        A_demand[cntry] = output.Aprime

    end

return Y, tradeflows, A_demand, hh, dist

end

# ##########################################################################
# ##########################################################################

function micro_trade_elasticity(R, W, p, home, source, model_params; tol_vfi = 1e-6, hh_solution_method = "nl-fixedpoint")

    hh = solve_household_problem(R, W, p, model_params, tol = tol_vfi, solution_method = hh_solution_method)

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

function compute_eq(R, W, p, model_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "nl-fixedpoint", stdist_sol_method = "nl-fixedpoint")
    # Does everything...
    # (1) Sovles hh problem
    # (2) Constructs stationary distribution

    hh = solve_household_problem(R, W, p, model_params, tol = tol_vfi, solution_method = hh_solution_method)

    dist = make_stationary_distribution(hh, model_params, tol = tol_dis, solution_method = stdist_sol_method)
    
return hh, dist

end

# ##########################################################################
# ##########################################################################

function make_stationary_distribution(household, model_params; tol = 1e-10, solution_method = "nl-fixedpoint") 

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

    λ = itterate_stationary_distribution(Q; tol = 1e-10)

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

function solve_household_problem(R, W, p, model_params; tol = 10^-6, solution_method = "nl-fixedpoint")

    if solution_method == "nl-fixedpoint"

        Kga, πprob, Tv = policy_function_fixedpoint(R, W, p, model_params; tol = tol)

    elseif solution_method == "itteration"

        Kga, πprob, Tv = policy_function_itteration(R, W, p, model_params; tol = tol, Niter = 500)
        
    end
    
    return household(Kga, πprob, Tv)

end

##########################################################################
##########################################################################

function policy_function_itteration(R, W, p, model_params; tol = 10^-6, Niter = 500)
    # this is the boiler plate vfi routine (1) make grid (2) itterate on 
    # bellman operator untill convergence. 
    #
    # as Fast/ ~faster than Matlab (but nothing is multithreaded here)
    # fastest is using nlsove fixed point to find situation where
    # v = bellman_operator(v)
    
    @unpack Na, Nshocks, Ncntry, statesize, β, σϵ = model_params

    gc = repeat(range(0.1,3,Na),1,Nshocks,Ncntry)
    v = -ones(Na, Nshocks, Ncntry)/(1-β)

    Kgc = similar(gc)
    Tv = similar(v)


    for iter in 1:Niter
        
        Kgc, Tv = coleman_operator(gc, v, R, W, p, model_params)[1:2]

        err = vec_max(Kgc, gc)

        if err < tol

            #println(err)

            break
        end

 
        copy!(gc, Kgc)

        copy!(v,Tv)

        if iter == Niter

          println("value function may not have converged")
          println("check the situation")
          
        end

    end

    Tv, Kga = coleman_operator(gc, v, R, W, p, model_params)[2:3]

    πprob = make_πprob(Tv, σϵ)

    return Kga, πprob, Tv
    
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

function policy_function_fixedpoint(R, W, p, model_params; tol = 10^-6)

    @unpack Na, Nshocks, Ncntry, statesize, β, σϵ = model_params

    #foo, foobar = policy_function_itteration(R, W, p, model_params, Niter = 2)
    #policy_o = vcat(foo, foobar)

    Kgc = repeat(range(0.1,3,Na),1,Nshocks,Ncntry)
    Tv = -ones(Na, Nshocks, Ncntry)/(1-β)

    policy_o = vcat(Kgc, Tv)

    # have seen some convergence issues sometimes

    K(policy) = coleman_operator(policy, R, W, p, model_params)

    solution = fixedpoint(K, policy_o, ftol = tol, method = :anderson);
    
    if solution.f_converged == false
        println("did not converge")
    end

    #Kgc = solution.zero[1:Na, :, :]

    #Tv = solution.zero[(Na+1):end, :, :]

    Tv, Kga = coleman_operator(solution.zero[1:Na, :, :], solution.zero[(Na+1):end, :, :], R, W, p, model_params)[2:3]

    πprob = make_πprob(Tv, σϵ)

    return Kga, πprob, Tv

end

##########################################################################
##########################################################################

function value_function_fixedpoint(R, W, p, model_params; tol = 10^-6)

    @unpack Ncntry, Na, Nshocks, β, mc, σϵ = model_params

    u = Array{eltype(R)}(undef, Na, Na, Nshocks, Ncntry)
    
    make_utility!(u, R, W, p, model_params)

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


##########################################################################
##########################################################################



