struct household{T}
    Tv::Array{T} # value function
    asset_policy::Array{T} # asset_policy
end

struct distribution{T}
    Q::Array{T} # transition matrix
    λ::Array{T} # λ
    state_index::Array{Tuple{Int64, Int64}} # index of states, lines up with λ
end

# ##########################################################################
# # functions used to find a solution to the households problem and then
# # the stationary equillibrium. 
# # Keep the idea in mind seperate model vs. solution technique.
# # so this files is about solution. Model enviornment is in envronment.jl file

# function clear_asset_market(Pces, W, τ_rev, R, model_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
#     vfi_solution_method = "nl-fixedpoint", stdist_sol_method = "nl-fixedpoint")

# @assert model_params.ϕ > 0.0

# hh = solve_household_problem(Pces, W, τ_rev, R, model_params; tol = tol_vfi, solution_method = vfi_solution_method)

# dist = make_stationary_distribution(hh, model_params, tol = tol_dis, solution_method = stdist_sol_method)

# output = aggregate(Pces, W, τ_rev, R, hh, dist, 1.0, model_params)

# return output.Aprime

# end

# function clear_asset_market(Pces, W, 
#     τ_rev, R::ForwardDiff.Dual, model_params; tol_vfi = 1e-6, tol_dis = 1e-10)
# # for some reason NLsolve fixedpoint does not work well with ForwardDiff.Dual numbers.
# # the work around is to use multiple dispatch. So if R a dual number...then revert to 
# # basic itterative methods to solve hh problem and stationary distribution.

# @assert model_params.ϕ > 0.0

# hh = solve_household_problem(Pces, W, τ_rev, R, model_params; tol = tol_vfi, solution_method = "vfi-itteration")

# dist = make_stationary_distribution(hh, model_params, tol = tol_dis, solution_method = "itteration")

# output = aggregate(Pces, W, τ_rev, R, hh, dist, 1.0, model_params)

# return output.Aprime

# end

# ##########################################################################
# ##########################################################################

# function compute_eq(Pces, W::ForwardDiff.Dual, τ_rev::ForwardDiff.Dual, 
#     R::ForwardDiff.Dual, model_params; tol_vfi = 1e-6, tol_dis = 1e-10)
# # Does everything...
# # (1) Sovles hh problem
# # (2) Constructs stationary distribution
# hh = solve_household_problem(Pces, W, τ_rev, R, model_params; tol = tol_vfi, solution_method = "vfi-itteration")

# dist = make_stationary_distribution(hh, model_params, tol = tol_dis, solution_method = "itteration")

# return hh, dist

# end

function compute_eq(W, R, model_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "nl-fixedpoint", stdist_sol_method = "nl-fixedpoint")
# Does everything...
# (1) Sovles hh problem
# (2) Constructs stationary distribution
    
    hh, Q, state_index = solve_household_problem(W, R, model_params, tol = tol_vfi, solution_method = hh_solution_method)

    dist = make_stationary_distribution(Q, state_index, model_params, tol = tol_dis, solution_method = stdist_sol_method)
    
return hh, dist

end

##########################################################################
##########################################################################

function solve_household_problem(W, R, model_params; tol = 10^-6, solution_method = "nl-fixedpoint")

    if solution_method == "nl-fixedpoint"

        solution = policy_function_fixedpoint(w, R, model_params; tol = tol)

        Kg = solution.zero

    elseif solution_method == "itteration"

        Kg = policy_function_itteration(w, R, model_params; tol = tol, Niter = 500)
        
    end
    
    @unpack Na, Nshocks, statesize = model_params
    
    state_index = Array{Tuple{eltype(Na), eltype(Na)}}(undef, statesize, 1)

    Q = Array{eltype(R)}(undef, statesize, statesize)
    
    Kg, asset_policy, Q, state_index, Tv = coleman_operator(Kg, Q, state_index, R, W, model_params);

    return household(Tv, asset_policy), Q, state_index

end

##########################################################################
##########################################################################

function make_stationary_distribution(Q, state_index, model_params; tol = 1e-10, solution_method = "nl-fixedpoint") 

@unpack Na, Nshocks, statesize = model_params

if solution_method == "nl-fixedpoint"

    g(λ) = law_of_motion(λ, transpose(Q))

    initialvalue = zeros(size(Q)[1], 1)

    initialvalue .= 1.0 / size(Q)[1]

    solution = fixedpoint(g, initialvalue, ftol = tol, method = :anderson)

    λ = solution.zero

elseif solution_method == "itteration"

    λ = itterate_stationary_distribution(Q; tol = 1e-10)

end

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
    
    for iter in 1:Niter
        
        #Lnew = transpose(Q) * λ
        
        Lnew = law_of_motion(λ, transpose(Q))
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

##########################################################################
##########################################################################

function policy_function_itteration(R, W, p, model_params; tol = 10^-6, Niter = 500)
    # this is the boiler plate vfi routine (1) make grid (2) itterate on 
    # bellman operator untill convergence. 
    #
    # as Fast/ ~faster than Matlab (but nothing is multithreaded here)
    # fastest is using nlsove fixed point to find situation where
    # v = bellman_operator(v)
    
    @unpack Na, Nshocks, Ncntry, statesize = model_params

    Q = Array{eltype(R)}(undef, statesize, statesize)

    gc = ones(Na, Nshocks, Ncntry)
    πprob = (1/Ncntry)*ones(Na, Nshocks, Ncntry)

    Kgc = similar(gc)
    Kπprob= similar(πprob)

    for iter in 1:Niter
        
        Kgc, Kπprob = coleman_operator(gc, πprob, Q, R, W, p, model_params)[1:2]
    
        err = maximum(abs, vcat(Kgc - gc, Kπprob - πprob) )

        err < tol && break
                
        gc = copy(Kgc)

        πprob = copy(Kπprob)

        if iter == Niter

          println("value function may not have converged")
          println("check the situation")
        end

    end

    Kgc, Kπprob, v, Q = coleman_operator(Kgc, Kπprob, Q, R, W, p, model_params)

    return Kgc, Kπprob, v, Q
    
end

##########################################################################
##########################################################################

function policy_function_fixedpoint(R, W, p, model_params; tol = 10^-6)

    @unpack Na, Nshocks, Ncntry, statesize = model_params

    Q = Array{eltype(R)}(undef, statesize, statesize)

    policyo = Array{Array{Float64}(undef, Na, Nshocks, Ncntry)}(undef,2)

    policyo[1] = ones(Na, Nshocks, Ncntry)
    policyo[2] = (1/Ncntry)*ones(Na, Nshocks, Ncntry)
  
    K(policy) = coleman_operator(policy, Q, R, W, p, model_params)

    solution = fixedpoint(K, policyo, ftol = tol, method = :anderson);
    
    if solution.f_converged == false
        println("did not converge")
    end

    Kgc, Kπprob, v, Q = coleman_operator(policy[1], policy[2], Q, R, W, p, model_params)

    return Kgc, Kπprob, v, Q

end

##########################################################################
##########################################################################

function value_function_fixedpoint(R, W, p, model_params; tol = 10^-6)

    @unpack Ncntry, Na, Nshocks, β, mc, σϵ = model_params

    u = Array{eltype(R)}(undef, Na, Na, Nshocks, Ncntry)
    
    make_utility!(u, W, R, p, model_params)

# define the inline function on the bellman operator. 
# so the input is v (other stuff is fixed)
    
    T(v) = bellman_operator_upwind(v, u, mc.p, β, σϵ)

    Vo = -ones(Na, Nshocks) / (1-β); #initial value
    Vo = convert(Array{eltype(u)}, Vo)

    solution = fixedpoint(T, Vo, ftol = tol, method = :anderson);
    
    if solution.f_converged == false
        println("did not converge")
    end

    return solution, u

end


##########################################################################
##########################################################################



