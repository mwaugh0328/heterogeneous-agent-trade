struct household{T}
    Tv::Array{T} # value function
    asset_policy::Array{T} # asset_policy
end

struct distribution{T}
    Q::Array{T} # transition matrix
    L::Array{T} # L
    state_index::Array{Tuple{Int64, Int64}} # index of states, lines up with L
end

##########################################################################
# functions used to find a solution to the households problem and then
# the stationary equillibrium. 
# Keep the idea in mind seperate model vs. solution technique.
# so this files is about solution. Model enviornment is in envronment.jl file

function clear_asset_market(Pces, W, τ_rev, R, model_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    vfi_solution_method = "nl-fixedpoint", stdist_sol_method = "nl-fixedpoint")

@assert model_params.ϕ > 0.0

hh = solve_household_problem(Pces, W, τ_rev, R, model_params; tol = tol_vfi, solution_method = vfi_solution_method)

dist = make_stationary_distribution(hh.asset_policy, model_params, tol = tol_dis, solution_method = stdist_sol_method)

aprime = get_aprime(hh.asset_policy, dist.state_index, model_params)

net_asset_demand = sum(aprime .* L, dims = 1)[1]

return net_asset_demand

end

##########################################################################
##########################################################################

function compute_eq(Pces, W, τ_rev, R, model_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    vfi_solution_method = "nl-fixedpoint", stdist_sol_method = "nl-fixedpoint")
# Does everything...
# (1) Sovles hh problem
# (2) Constructs stationary distribution
hh = solve_household_problem(Pces, W, τ_rev, R, model_params; tol = tol_vfi, solution_method = vfi_solution_method)

dist = make_stationary_distribution(hh.asset_policy, model_params, tol = tol_dis, solution_method = stdist_sol_method)

return hh, dist

end

##########################################################################
##########################################################################

function solve_household_problem(Pces, W, τ_rev, R, model_params; tol = 10^-6, solution_method = "nl-fixedpoint")

    if solution_method == "nl-fixedpoint"

        solution, u = value_function_fixedpoint(Pces, W, τ_rev, R, model_params; tol = tol)

        Tv = solution.zero

    elseif solution_method == "vfi-itteration"

        Tv, u = value_function_itteration(Pces, W, τ_rev, R, model_params; tol = tol, Niter = 500)
        
    end
    
    @unpack Na, Nshocks,
     mc, β, σa = model_params

     Tv, asset_policy = bellman_operator_policy(Tv, u, mc.p, β, σa);

    return household(Tv, asset_policy )

end

##########################################################################
##########################################################################

function make_stationary_distribution(asset_policy, model_params; tol = 1e-10, solution_method = "nl-fixedpoint") 

@unpack Na, Nshocks, statesize = model_params

state_index = Array{Tuple{eltype(Na), eltype(Na)}}(undef, statesize, 1)

Q = Array{eltype(asset_policy)}(undef, statesize, statesize)

make_Q!(Q, state_index, asset_policy, model_params)

if solution_method == "nl-fixedpoint"

    g(L) = law_of_motion(L, transpose(Q))

    initialvalue = zeros(size(Q)[1], 1)

    initialvalue .= 1.0 / size(Q)[1]

    solution = fixedpoint(g, initialvalue, ftol = tol, method = :anderson)

    stationary_distribution = solution.zero

elseif solution_method == "itteration"

    stationary_distribution = itterate_stationary_distribution(Q; tol = 1e-10)

end

return distribution(Q, stationary_distribution, state_index)

end

##########################################################################
##########################################################################
function itterate_stationary_distribution(Q; tol = 1e-10, Niter = 5000) 
    # this is faster than the quant econ canned routine
    # from lyon-waugh implementation. Not as fast as using
    # NLsolve fixedpoint routine below.
    
    # Takes the transition matrix above then we know a stationary distribution
    # must satisfy the fixed point relationsip L = Q'*L  
    
    L = zeros(size(Q)[1], 1)
    
    L[1] = 1.0
    
    Lnew = similar(L)
    
    for iter in 1:Niter
        
        #Lnew = transpose(Q) * L
        
        Lnew = law_of_motion(L, transpose(Q))
        #this ordering is also better in julia
        # than my matlab implementation of Q*L (1, na*nshock)
                
        err = maximum(abs, L - Lnew)
        
        copy!(L, Lnew)
        # this surprisingly makes a big difference
        # but in the vfi it causes a slowdown?
        
        err < tol && break

        if iter == Niter

            println("distribution may not have converged")
            println("check the situation")
    
        end
        
    end

    return L

end

##########################################################################
##########################################################################


function value_function_itteration(Pces, W, τ_rev, R, model_params; tol = 10^-6, Niter = 500)
    # this is the boiler plate vfi routine (1) make grid (2) itterate on 
    # bellman operator untill convergence. Policy 
    # functions are then outputed as cartesian index
    #
    # as Fast/ ~faster than Matlab (but nothing is multithreaded here)
    # fastest is using nlsove fixed point to find situation where
    # v = bellman_operator(v)
    
    @unpack Na, Nshocks, β, mc = model_params

    u = Array{eltype(R)}(undef, Na, Na, Nshocks)

    make_utility!(u, Pces, W, τ_rev, R, model_params)

    v = -ones(Na, Nshocks) / (1-β); #initial value
    
    v = convert(Array{eltype(u)}, v)
    Tv = similar(v)

    for iter in 1:Niter
        
        Tv = bellman_operator_upwind(v, u, mc.p, β) 
        #there is some advantage of having it
        # explicity, not always recreating the Tv 
        # array in the function
    
        err = maximum(abs, Tv - v)

        err < tol && break
                
        v = copy(Tv)

        if iter == Niter

          println("value function may not have converged")
          println("check the situation")
        end

    end

    return Tv, u
    
end

##########################################################################
##########################################################################

function value_function_fixedpoint(Pces, W, τ_rev, R, model_params; tol = 10^-6)

    @unpack Na, Nshocks, β, mc = model_params

    u = Array{eltype(R)}(undef, Na, Na, Nshocks)

    make_utility!(u, Pces, W, τ_rev, R, model_params)

# define the inline function on the bellman operator. 
# so the input is v (other stuff is fixed)
    
    T(v) = bellman_operator_upwind(v, u, mc.p, β)

    Vo = -ones(Na, Nshocks) / (1-β); #initial value

    solution = fixedpoint(T, Vo, ftol = tol, method = :anderson);
    
    if solution.f_converged == false
        println("did not converge")
    end

    return solution, u

end


##########################################################################
##########################################################################

