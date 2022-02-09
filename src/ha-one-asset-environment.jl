using BenchmarkTools, SpecialFunctions
#LinearAlgebra.BLAS.set_num_threads(12)
using Parameters
using QuantEcon
using BasisMatrices
using NLsolve
using Octavian

##########################################################################

model_params = @with_kw (
    β= 0.96,
    γ = 1.0,
    ϕ = 1.0,
    amax = 3.0,
    Na= 50,
    agrid = range(-ϕ, amax, length = Na),
    Nshocks = 5,
    statesize = Int(Na*Nshocks),
    ρ = 0.90,
    σ = 0.10,
    mc = tauchen(Nshocks, ρ, σ),
    σa = 0.01,
)

##########################################################################
function bellman_operator_policy(v, u, mc, β, σa) 
    # the value function /bellman operator that takes a v then
    # returns a Tv and policy functions.
    # The asset policy is smoothed and are probabilities over a' tomorrow. 

    Na = size(u)[1]
    Nshocks = size(u)[3]

    asset_policy = Array{eltype(u)}(undef, Na, Na,  Nshocks)

    Tv = Array{eltype(u)}(undef, Na, Nshocks)

    for shockstate = 1:Nshocks

        βEV = compute_EV(β*v, mc[shockstate, :])

        foo = u[:, :, shockstate] .+ βEV

        maximum!(view(Tv, : , shockstate), foo )

        foo .-= ( Tv[ :, shockstate]) 
        # this operation is creating NaN's in foo
        # then propogates.  max is -Inf and then 
        # the subtraction of -Inf so (-Inf - (-Inf) = NaN)
                    
        foo[isnan.(foo)] .= 0.0
                    
        @views asset_policy[ :, :, shockstate] = exp.( foo ./ σa ) ./ sum( exp.( foo ./ σa ) , dims = 2) 
                   
    end
 
    return Tv, asset_policy
    
end

##########################################################################
##########################################################################

function bellman_operator(v, u, mc, β)
    # basic value function /bellman operator that takes a v then
    # returns a TV. 
    # v and Tv are setup each individual entry has v(a,z) as in 
    # model.pdf notes.
           
    Na = size(u)[1]
    Nshocks = size(u)[3]

    # stores the value function Tv_{j}(a,z)
    Tv = Array{eltype(u)}(undef, Na, Nshocks)

    # eltype(u) ensures that Tvmove is same type as u

    for shockstate = 1:Nshocks
        # work through each shock state

        βEV = compute_EV(β*v, mc[shockstate, :])
             # Compute expected value 
                
        maximum!(view(Tv, : , shockstate), u[:, :, shockstate] .+ βEV )

    end
    # this returns a array of Na, Nschocks, Ncars, 1 
    # so need drop the last 1 dimension...
    
    return Tv
    
end



##########################################################################
##########################################################################

function bellman_operator_upwind(v, u, mc, β) 
    # the value function /bellman operator that takes a v then
    # returns a TV.
    # 
    # upwind approach... so the EV is evaluated with TV, not V...
    # usually yeilds faster convergence
    #
    # v and Tv are setup each individual entry has v(a,z) as in 
    # model.pdf notes.

    Nshocks = size(u)[3]
    
    Tv = copy(v)

    for shockstate = 1:Nshocks
        # work through each shock state

        @views βEV = compute_EV(β*Tv, mc[shockstate, :])
        # Compute expected value, key difference here is
        # Tv is used to compute EV, not gueessed v.  
                
        maximum!(view(Tv, : , shockstate), u[:, :, shockstate] .+ βEV )

    end
    
    return Tv
      
end


##########################################################################
##########################################################################
function compute_EV(v, mc_probs)
    # construct expected value
    # given the transition probabiltes from the markov chain...
    
    # so mc_probs should be a 1 by Nshock row vector 
    # then v should be Nassset state, Nshock state matrix
    # we want to integrate this so transpose, so for a given
    # asset state integrate accross different shock outcomes
    
    #return  (mc_probs' * v' )

    return matmul(mc_probs', v' )
    
    # updated: for some reason matmul in this place does not work with autodiff
    # not sure why, it works with autodiff in the construction of stationary distribution 
    
    # uses Octavian a pure julia implementation of matrix multiplication
    # standard * calls BLAS C(?) wrapped implementation
    # seems faster than *
end

##########################################################################
##########################################################################
function make_utility!(utility_grid, Pces, W, τ_rev, R, model_params) 
    # take prices and model parameters and returns utility function
    # R is gross real interest rate ∈ (β, 1/β)
    # W is wage per effeciency unit
    
    @unpack Na, Nshocks, mc, agrid, γ = model_params
    
    a =  reshape(agrid, Na, 1)
    #assets today
    
    a_prime = transpose(a)
    #assets tomorrow

    shock_level = exp.(mc.state_values)

    for shockstate = 1:Nshocks
    #shock state

        wz = labor_income( shock_level[shockstate] , W )

        c = consumption(Pces, τ_rev, R*a, a_prime, wz)
                # takes assets states, shock state -> consumption from
                # budget constraint
        
        utility_grid[:, :, shockstate] = utility.(c, γ)
        
    end
    
end


##########################################################################
##########################################################################
function labor_income(shock, W)
    # computes labor income. it's simple now, but need so it can be 
    # more complicated later
    # W is wage per effeciency unit

    return shock * W

end
            
##########################################################################
##########################################################################

function consumption(Pces, τ_rev, Ra, ap, wz)
        # constructs consumption via the hh budget constraint
        # .- boadcasting is used to turn vectors (a,ap) into a grid
        # over c (na by na)
    
        return @. Ra - ap + (wz / Pces)  + (τ_rev)
        
end
    
##########################################################################
##########################################################################
function log_utility(c)

    (c < 1e-10 ? -Inf : @fastmath log(c) )

end

##########################################################################
##########################################################################

function utility(c, γ)
    # maps consumption into utility with the CRRA specification
    # log it γ is close to one
    
    if γ == 1.0
        
        (c < 1e-10 ? -Inf : @fastmath log(c) )

    else
        (c < 1e-10 ? -Inf : @fastmath c^( 1.0 - γ) / (1.0 - γ))

    end

end

##########################################################################
##########################################################################

function law_of_motion(L , Q_tran)
    # Takes the transition matrix Q and given a distribution L, advances it
    # where Lnew = Q'*L
    # usefull to construct stationary where L = Q'*L
    
    #return Q_tran * L
    
    return matmul(Q_tran, L) 
    # this is using Octavian a fast, pure julia package for matrix multiplicaiton
    # if LinearAlgebra.BLAS.set_num_threads(Threads.nthreads()) need to be carefull
    # as it is chosing # threads independently. When free, matmul beats by half, with 24 
    # threads its close. With 1 it's 1/4 faster.

    # note that this works with autodiff here
    
end

##########################################################################
##########################################################################

function make_Q!(Q, state_index, asset_policy, model_params)

    @unpack Na, Nshocks, mc= model_params
    
    fill!(Q, 0.0) # this is all setup assumeing Q is zero everywehre

    @inbounds for shk = 1:Nshocks
    
        shk_counter = Int((shk - 1)*Na)

        for ast = 1:Na

            today = ast + shk_counter 

            state_index[today] = (ast, shk)

            for shkprime = 1:Nshocks
    
                shk_counter_prime = Int((shkprime - 1)*Na)

                for astprime = 1:Na

                    tommorow = astprime + shk_counter_prime 

                    Q[today, tommorow] += mc.p[shk, shkprime] * asset_policy[ast, astprime, shk]

                end

            end

        end

    end

end
