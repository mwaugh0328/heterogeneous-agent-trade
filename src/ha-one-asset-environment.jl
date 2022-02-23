using BenchmarkTools, SpecialFunctions
#LinearAlgebra.BLAS.set_num_threads(12)
using Parameters
using QuantEcon
using BasisMatrices
using NLsolve
using Octavian
using ForwardDiff

##########################################################################

@with_kw struct model_params
    β::Float64 = 0.96
    γ::Float64 = 1.0
    ϕ::Float64 = 1.0
    amax::Float64 = 5.0
    Na::Int64 = 100
    agrid::Array{Float64, 1} = convert(Array{Float64, 1}, range(-ϕ, amax, length = Na))
    Nshocks::Int64 = 5
    statesize::Int64 = Int(Na*Nshocks)
    ρ::Float64 = 0.90
    σ::Float64 = 0.10
    mc::MarkovChain{Float64, Matrix{Float64}, 
    StepRangeLen{Float64, Base.TwicePrecision{Float64},
     Base.TwicePrecision{Float64}, Int64}} = tauchen(Nshocks, ρ, σ)
    σa::Float64 = 0.01
    σw::Float64 = 0.25
    ϑ::Float64 = 1.5
    Woptions::Int64 = 2
end



##########################################################################
function bellman_operator_policy(v, u, mc, β, σa, σw) 
    # the value function /bellman operator that takes a v then
    # returns a Tv and policy functions.
    # The asset policy is smoothed and are probabilities over a' tomorrow. 

    Na = size(u)[1]
    Nshocks = size(u)[3]
    Woptions = size(u)[4]

    asset_policy = Array{eltype(u)}(undef, Na, Na,  Nshocks, Woptions)
    work_policy = Array{eltype(u)}(undef, Na, Nshocks, Woptions)
    foo = Array{eltype(u)}(undef, Na, Na)

    Tvwork = Array{eltype(u)}(undef, Na, Nshocks, Woptions)
    Tv = Array{eltype(u)}(undef, Na, Nshocks)

    @inbounds @views for shockstate = 1:Nshocks

        βEV = compute_EV(β*v, mc[shockstate, :])

        for wrkchoice = 1:Woptions

            foo .= u[:, :, shockstate, wrkchoice] .+ βEV

            maximum!(view(Tvwork, : , shockstate, wrkchoice), foo )

            foo .-= ( Tvwork[ :, shockstate, wrkchoice]) 
            foo[isnan.(foo)] .= 0.0
            # this operation is creating NaN's in foo
            # then propogates.  max is -Inf and then 
            # the subtraction of -Inf so (-Inf - (-Inf) = NaN)

            asset_policy[ :, :, shockstate, wrkchoice] = exp.( foo ./ σa ) ./ sum( exp.( foo ./ σa ) , dims = 2) 

        end

        make_Tv!(view(Tv, :, shockstate), view(Tvwork,:,shockstate,:), view(work_policy,:,shockstate,:), σw)
               
    end

    return household(Tv, asset_policy, work_policy )
    
end

##########################################################################
##########################################################################

function bellman_operator(v, u, mc, β, σw)
    # basic value function /bellman operator that takes a v then
    # returns a TV. 
    # v and Tv are setup each individual entry has v(a,z) as in 
    # model.pdf notes.
           
    Na = size(u)[1]
    Nshocks = size(u)[3]
    Woptions = size(u)[4]

    # stores the value function Tv_{j}(a,z)
    
    Tv = Array{eltype(u)}(undef, Na, Nshocks)
    Tvwork = Array{eltype(u)}(undef, Na, Nshocks, Woptions)

    @inbounds @views for shockstate = 1:Nshocks
        # work through each shock state

        # Compute expected value 
        βEV = compute_EV(β*v, view(mc,shockstate, :))

        for wrkchoice = 1:Woptions

            # find best value for each choice
            maximum!(view(Tvwork, : , shockstate, wrkchoice), u[:, :, shockstate, wrkchoice] .+ βEV )

        end

        make_Tv!(view(Tv,:,shockstate), view(Tvwork, :, shockstate,:), σw)

    end

    return Tv
    
end




##########################################################################
##########################################################################

function bellman_operator_upwind(v, u, mc, β, σw) 
    # the value function /bellman operator that takes a v then
    # returns a TV.
    # 
    # upwind approach... so the EV is evaluated with TV, not V...
    # usually yeilds faster convergence
    #
    # v and Tv are setup each individual entry has v(a,z) as in 
    # model.pdf notes.
    Na = size(u)[1]
    Nshocks = size(u)[3]
    Woptions = size(u)[4]
    
    Tv = copy(v)

    Tvwork = Array{eltype(u)}(undef, Na, Nshocks, Woptions)

    @inbounds @views for shockstate = 1:Nshocks
        # work through each shock state

        # Compute expected value 
         βEV = compute_EV(β*Tv, mc[shockstate, :])

        for wrkchoice = 1:Woptions

            # find best value for each choice
            maximum!(view(Tvwork, : , shockstate, wrkchoice), u[:, :, shockstate, wrkchoice] .+ βEV )

        end

        make_Tv!(view(Tv,:,shockstate), view(Tvwork, :, shockstate,:), σw)

    end

    return Tv
      
end

##########################################################################
##########################################################################

function make_Tv!(Tv, Tvwork, σw)
    #makes the Tv function with the scaling by max across options
    # this is abit faster than old version as it eliminates need to
    # create new variables...

    Tvworkmax = maximum(Tvwork, dims = 2)

    Tvwork .-= Tvworkmax # broadcast, inplace subtract off max value

    Tvwork[isnan.(Tvwork)] .= 0.0

    Tv .= σw.*log.( sum( exp.(  Tvwork ./ σw ) , dims = 2) )  .+ Tvworkmax

end

function make_Tv!(Tv, Tvwork, work_policy, σw)
    #multiple dispact version to fill in policy function
    #makes the Tv function with the scaling by max across options
    # this is abit faster than old version as it eliminates need to
    # create new variables...

    Tvworkmax = maximum(Tvwork, dims = 2)

    Tvwork .-= Tvworkmax # broadcast, inplace subtract off max value

    Tvwork[isnan.(Tvwork)] .= 0.0

    work_policy .= exp.( Tvwork ./ σw ) ./ sum( exp.( Tvwork ./ σw ) , dims = 2) 

    Tv .= σw.*log.( sum( exp.(  Tvwork ./ σw ) , dims = 2) )  .+ Tvworkmax

end


##########################################################################
##########################################################################
function compute_EV(v::Array{Float64}, mc_probs::Array{Float64})
    # construct expected value
    # given the transition probabiltes from the markov chain...
    
    # so mc_probs should be a 1 by Nshock row vector 
    # then v should be Nassset state, Nshock state matrix
    # we want to integrate this so transpose, so for a given
    # asset state integrate accross different shock outcomes
    
    return matmul(mc_probs', v' )

    #matmul is from Octavian packaage, pure julia matrix multiplicaiton
    # can be meaningfully faster than BLAS
    # problem is it does not play well with ForwardDiff when mc_probs
    # is a subarray. So solution is to use mulitple dispactch. IF anything,
    # then use BLAS. If Array then specilize and do matmul

end

##########################################################################
##########################################################################
function compute_EV(v, mc_probs) 
    # multiple dispatch version if v is a dual number

    return  (mc_probs' * v' )
        
    # updated: for some reason matmul in this place does not work with autodiff
    # not sure why, it works with autodiff in the construction of stationary distribution 

end

##########################################################################
##########################################################################
function make_utility!(utility_grid, Pces, W, τ_rev, R, model_params) 
    # take prices and model parameters and returns utility function
    # R is gross real interest rate ∈ (β, 1/β)
    # W is wage per effeciency unit
    
    @unpack Na, Nshocks, Woptions, mc, agrid, γ, ϑ = model_params
    
    a =  reshape(agrid, Na, 1)
    #assets today
    
    a_prime = transpose(a)
    #assets tomorrow

    shock_level = exp.(mc.state_values)

    @inbounds @views for shockstate = 1:Nshocks
    #shock state

        for workchoice = 1:Woptions
        
            wz = labor_income( shock_level[shockstate] , W, workchoice)

            c = consumption(Pces, τ_rev, R.*a, a_prime, wz)
                # takes assets states, shock state -> consumption from
                # budget constraint
        
            utility_grid[:, :, shockstate, workchoice] .= utility.(c, γ, ϑ, workchoice)

        end
        
    end
    
end


##########################################################################
##########################################################################
function labor_income(shock, W, workchoice)
    # computes labor income. it's simple now, but need so it can be 
    # more complicated later
    # W is wage per effeciency unit

    if workchoice == 1
        #you are working

        return shock * W

    else
        # not working get nothing

        return 0.0

    end

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
function utility(c, γ)
    # maps consumption into utility with the CRRA specification
    # log it γ is close to one

    if γ == 1.0
        
        (c < 1e-10 ? -Inf : log(c) )

    else
        (c < 1e-10 ? -Inf : c^( 1.0 - γ) / (1.0 - γ))

    end

end

# Using mulitiple dispatch here...

function utility(c, γ, ϑ, workchoice)

    if workchoice == 1

        utility(c, γ) # don't enjoy it

    else

        utility(c, γ) + ϑ # enjoy leisure

    end

end

##########################################################################
##########################################################################

function law_of_motion(λ , Q_tran)
    # Takes the transition matrix Q and given a distribution λ, advances it
    # where Lnew = Q'*λ
    # usefull to construct stationary where λ = Q'*λ
    
    #return Q_tran * λ
    
    return matmul(Q_tran, λ) 
    # this is using Octavian a fast, pure julia package for matrix multiplicaiton
    # if LinearAlgebra.BLAS.set_num_threads(Threads.nthreads()) need to be carefull
    # as it is chosing # threads independently. When free, matmul beats by half, with 24 
    # threads its close. With 1 it's 1/4 faster.

    # note that this works with autodiff here
    
end

##########################################################################
##########################################################################

function make_Q!(Q, state_index, asset_policy, work_policy, model_params)

    @unpack Na, Nshocks, Woptions, mc= model_params
    
    fill!(Q, 0.0) # this is all setup assumeing Q is zero everywehre

    @inbounds @views for shk = 1:Nshocks
    
        shk_counter = Int((shk - 1)*Na)

        for ast = 1:Na

            today = ast + shk_counter 

            state_index[today] = (ast, shk)

            for shkprime = 1:Nshocks
    
                shk_counter_prime = Int((shkprime - 1)*Na)

                for astprime = 1:Na

                    tommorow = astprime + shk_counter_prime 

                    for wrk = 1:Woptions

                        Q[today, tommorow] += mc.p[shk, shkprime] * asset_policy[ast, astprime, shk, wrk] * work_policy[ast, shk, wrk]

                    end

                end

            end

        end

    end

end
