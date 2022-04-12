using MKL # Bug in Julia and openBlas, #Intel version
using BenchmarkTools, SpecialFunctions
#LinearAlgebra.BLAS.set_num_threads(12)
using Parameters
using QuantEcon
using BasisMatrices
using NLsolve
using Octavian
using ForwardDiff
using Interpolations
using LinearAlgebra


##########################################################################

@with_kw struct model_params
    β::Float64 = 0.96
    γ::Float64 = 2.0
    ϕ::Float64 = 0.0
    amax::Float64 = 14.0
    Ncntry::Int64 = 2
    Na::Int64 = 50
    agrid::Array{Float64, 1} = convert(Array{Float64, 1}, range(-ϕ, amax, length = Na))
    Nshocks::Int64 = 2
    statesize::Int64 = Int(Na*Nshocks)
    ρ::Float64 = 0.60
    σ::Float64 = 0.30
    mc::MarkovChain{Float64, Matrix{Float64}, 
    StepRangeLen{Float64, Base.TwicePrecision{Float64},
     Base.TwicePrecision{Float64}, Int64}} = tauchen(Nshocks, ρ, σ)
end

##########################################################################
##########################################################################



function coleman_operator_test(c, πprob, R, W, p, model_params)
    # policyfun = 

    @unpack agrid, mc, β, γ, Nshocks, Ncntry = model_params

    shocks = exp.(mc.state_values)

    Kg = similar(c)

    muc_ϵ = Array{Float64}(undef, Na, Nshocks)
    
    #first step is to integrate out future taste shock

    ∑π_ϵ!(muc_ϵ, c, πprob, p, γ)
    # this is the ∑ π_j(a,z) * muc_j(a,z) / p_j

    Emuc = β*R*( matmul( muc_ϵ , mc.p'))
    # this integrates over z

    for cntry = 1:Ncntry

        gc = muc_inverse.( p[cntry] * Emuc, γ)
        # invert from rhs of euler equation
        # key is how price shows up here.
        
        ã = (p[cntry] * gc .+ agrid .- W*shocks') / R
        # off budget constraint a = (p_jc_j + a' - w*z ) / R

        for nshk = 1:Nshocks
    
            foo = LinearInterpolation(ã[:,nshk], agrid, extrapolation_bc = (Flat(), Flat()) )
    
            Kg[:, nshk, cntry] = ( -foo.(agrid) .+ R*agrid .+ W*shocks[nshk] ) / p[cntry]
            # again off budget constraint pc = -a + Ra + wz
    
        end

    end

    return Kg

end


function ∑π_ϵ!(muc_ϵ, c, πprob, p, γ)

    Na, Nshocks = size(c)[1:2]

    for ast = 1:Na

        for shk = 1:Nshocks
            
            muc_ϵ[ast,shk] =  dot((πprob[ast, shk, :]), ( muc.(c[ast, shk, :] , γ) ./ p ))
            # this is the inside part of 54
            # not sure why matmul does not work here...kicks back error

        end
    end

end

##########################################################################
##########################################################################

function coleman_operator(c, R, W, model_params)
    # mulitple-disptach here, when c is passed it returns Kgc

    @unpack agrid, mc, β, γ, Nshocks = model_params

    shocks = exp.(mc.state_values)

    Kg = similar(c)

    gc = muc_inverse.( β*R*( matmul( muc.(c, γ) , mc.p')), γ) 

    ã = (gc .+ agrid .- W*shocks') / R

    for nshk = 1:Nshocks

        foo = LinearInterpolation(ã[:,nshk], agrid, extrapolation_bc = (Flat(), Linear()) )

        Kg[:, nshk] = -foo.(agrid) .+ R*agrid .+ W*shocks[nshk]

    end

    return Kg

end

##########################################################################
##########################################################################

function coleman_operator(c, Q, state_index, R, W, model_params)
    # multiple dispatch here, returns Kg, gq, and v

    @unpack agrid, mc, β, γ, Na, Nshocks, statesize = model_params

    shocks = exp.(mc.state_values)

    v = similar(c)
    aprime = similar(c)
    u = similar(c)
    aindex = Array{Int64}(undef, size(c))

    Kg = similar(c)

    # Step 1: Infer consumption today, given policy function, tommorow

    gc = muc_inverse.( β*R*( muc.(c, γ) * mc.p'), γ) 
    # muc = β*R*E(muc') 
    # and then to get consumption we have
    # c(ã,z) = muc_inverse ( β*R*E(muc(a',z')) ) 
    # so this tells us consumption today, given that we had state
    # a' tommorow (which is on the grid)

    # Step 2: Infer the state associated with that consumption today.

    ã = (gc .+ agrid .- W*shocks') / R

    # the budget constraint is
    # c + a' = R*a + w*z
    # so what we have inferd is c(ã,z) and an associated a',
    # but we don't know what ã is, so back it out from the budget constraint.
    # and now we have a map from states ã,z -> consumption and a'

    # Step 3: Push it back on to the grid...

    for nshk = 1:Nshocks

        foo = LinearInterpolation(ã[:,nshk], agrid, extrapolation_bc = (Flat(), Flat()) )
        # must have extrapolation here, otherwise an error is thrown when extrapolation 
        # is needed, Flat forces it to be at the bound. Mechanichally, what this does is that
        # if outside the bounds (defined by the grid), it extrapolates to the bound. 
        # So if ã is below the borrowing constraint -> the extrapolated value is ϕ

        # this creates a function "foo" which is created by interpoliating on 
        # a map from ã into a' where a' is on the grid.
        # then figure out what c(a,z) is by evaluating this map at a where a
        # is on the grid

        aprime[:, nshk] = foo.(agrid)

        Kg[:, nshk] = -aprime[:, nshk] .+ R*agrid .+ W*shocks[nshk]
        # this is the last step.
        # c = -a' + R*a + w*z

        u[:, nshk] = utility.(Kg[:, nshk], γ)

        for ast = 1:Na

            aindex[ast,nshk] = find_asset_position(aprime[ast, nshk], agrid)

        end

    end

    make_Q!(Q, state_index, aprime, aindex, model_params)

    v = (lu(I - β*Q)) \ u[:]
    # question is what is fastest way? # Open BLAS issues when big
    # a discourse happend to suggest lu above
    
    v = reshape(v, Na, Nshocks)

    return Kg, aprime, Q, state_index, v

end

##########################################################################
##########################################################################

function find_asset_position(ga, agrid)

    return findfirst(x -> (x >= ga), agrid)

end

##########################################################################
##########################################################################

function make_Q!(Q, state_index, asset_policy, asset_policy_index, model_params)
    # this is not optimized to exploit column major order.

    @unpack Na, Nshocks, mc, agrid = model_params
    
    fill!(Q, 0.0) # this is all setup assumeing Q is zero everywehre

    @inbounds @views for shk = 1:Nshocks
    
        shk_counter = Int((shk - 1)*Na)

        for ast = 1:Na

            today = ast + shk_counter 

            state_index[today] = (ast, shk)

            aprime_h = asset_policy_index[ast, shk]
            aprime_l = max(asset_policy_index[ast, shk] - 1, 1)
            
            p = 1.0 - (asset_policy[ast,shk] - agrid[aprime_l]) / (agrid[aprime_h] - agrid[aprime_l])

            if isnan(p)
                p = 1.0
            end

            for shkprime = 1:Nshocks
    
                shk_counter_prime = Int((shkprime - 1)*Na)

                for astprime = 1:Na

                    tommorow = astprime + shk_counter_prime 

                    if astprime == aprime_l

                        Q[today, tommorow] = ( p )*mc.p[shk, shkprime]

                    elseif astprime == aprime_h

                        Q[today, tommorow] = ( 1.0 - p )*mc.p[shk, shkprime]

                    end

                end

            end

        end

    end

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
    
    return matmul(mc_probs', v' )

    #matmul is from Octavian packaage, pure julia matrix multiplicaiton
    # can be meaningfully faster than BLAS
    # problem is it does not play well with ForwardDiff when mc_probs
    # is a subarray. So solution is to use mulitple dispactch. IF anything,
    # then use BLAS. If Array then specilize and do matmul

end

##########################################################################
##########################################################################

function muc(c, γ)
    # the marginal utility of consumption with
    # CRRA preferences

    c^(-γ) 

end

##########################################################################
##########################################################################

function muc_inverse(c, γ)
     # the marginal utility of consumption with
    # CRRA preferences

    c^( 1.0 / -γ )

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

##########################################################################
##########################################################################

function make_utility!(utility_grid, W, R, model_params) 
    # take prices and model parameters and returns utility function
    # R is gross real interest rate ∈ (β, 1/β)
    # W is wage per effeciency unit
    
    @unpack Na, Nshocks, mc, agrid, γ = model_params
    
    a = copy(agrid)
    #assets today
    
    a_prime = transpose(a)
    #assets tomorrow

    shock_level = exp.(mc.state_values)

    @inbounds @views for shockstate = 1:Nshocks
    #shock state

        wz = labor_income( shock_level[shockstate] , W)

        c = consumption(R.*a, a_prime, wz)
                # takes assets states, shock state -> consumption from
                # budget constraint
        
        utility_grid[:, :, shockstate] .= utility.(c, γ)
        
    end
    
end

##########################################################################
##########################################################################

function consumption(Ra, ap, wz)
    # constructs consumption via the hh budget constraint
    # .- boadcasting is used to turn vectors (a,ap) into a grid
    # over c (na by na)

    return @. Ra - ap + wz  
    
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

function bellman_operator_upwind(v, u, mc, β) 
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

    Tv = copy(v)

    for shockstate = 1:Nshocks
        # work through each shock state

        # Compute expected value 
         βEV = compute_EV(β*Tv, mc[shockstate, :])

        maximum!(view(Tv, : , shockstate), u[:, :, shockstate] .+ βEV )

    end

    return Tv
      
end