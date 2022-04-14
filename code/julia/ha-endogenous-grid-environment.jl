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
    σϵ::Float64 = 0.25
    Na::Int64 = 50
    agrid::Array{Float64, 1} = convert(Array{Float64, 1}, range(-ϕ, amax, length = Na))
    Nshocks::Int64 = 2
    statesize::Int64 = Int(Na*Nshocks*Ncntry)
    ρ::Float64 = 0.60
    σ::Float64 = 0.30
    mc::MarkovChain{Float64, Matrix{Float64}, 
    StepRangeLen{Float64, Base.TwicePrecision{Float64},
     Base.TwicePrecision{Float64}, Int64}} = tauchen(Nshocks, ρ, σ)
end

##########################################################################
##########################################################################

function coleman_operator(c, πprob, Q, R, W, p, model_params)
    # todo, better setup input of policy which is 
    # a consumption function and then choice probability

    ######################################################################
    # Organization 

    @unpack agrid, mc, β, γ, σϵ, Na, Nshocks, Ncntry, statesize = model_params

    shocks = exp.(mc.state_values)

    u = similar(c)
    π_selection = similar(c)
    aprime = similar(c)
    v = similar(c)
    aindex = Array{Int64}(undef, size(c))

    Kg = similar(c)
    Kπprob = similar(πprob)

    muc_ϵ = Array{eltype(R)}(undef, Na, Nshocks)

    ######################################################################
    # Now implement EGM method...

    #Step (1) is to integrate out future taste shock

    ∑π_ϵ!(muc_ϵ, c, πprob, p, γ)
    # this is the ∑ π_j(a,z) * muc_j(a,z) / p_j

    #Step (2) is to integrate out z'

    Emuc = β*R*( matmul( muc_ϵ , mc.p'))
    # this integrates over z

    #Step (3) Work through each county option
    @inbounds @views for cntry = 1:Ncntry

        gc = muc_inverse.( p[cntry] * Emuc, γ)
        # invert from rhs of euler equation
        
        ã = (p[cntry] * gc .+ agrid .- W*shocks') / R
        # off budget constraint a = (p_jc_j + a' - w*z ) / R

        # then linear interpolation to get back on grid.
        for shk = 1:Nshocks
    
            foo = LinearInterpolation(ã[:,shk], agrid, extrapolation_bc = (Flat(), Flat()) )

            aprime[:, shk, cntry] = foo.(agrid)
    
            Kg[:, shk, cntry] = ( -aprime[:, shk, cntry] .+ R*agrid .+ W*shocks[shk] ) / p[cntry]
            # again off budget constraint pc = -a + Ra + wz

            π_selection[:, shk, cntry] = -σϵ*log.(πprob[:, shk, cntry])

            u[:, shk, cntry] = utility.(Kg[:, shk, cntry], γ) + π_selection[:, shk, cntry]
            # the last bit factors in expected value of preference shock

            for ast = 1:Na

                aindex[ast, shk, cntry] = find_asset_position(aprime[ast, shk, cntry], agrid)

            end
    
        end

    end

    make_Q!(Q, aprime, aindex, πprob, model_params)
    # Think of the state as (a,z,j) so Q makes the transition matrix from
    # (a,z,j) -> (a',z',j') 

    vj = (lu(I - β*Q)) \ vec(u)
    # then given Q and note how u is defined over (a,z,j) we just apply the 
    # same procedure giving v(a,z,j) in each position
    #vj = itterate_v(vec(u),  β, Q)

    vj = reshape(vj , Na, Nshocks, Ncntry)
    # reshape v

    Kπprob = make_πprob(vj .- π_selection, σϵ)
    # this computes updated choice probablities given v
    # Important! need to net off selection correction

    v = reshape( sum( vj.*Kπprob, dims = 3), Na, Nshocks)
    # then this computes the expected v prior to prefernce shock
    # realization. not log sum thing because E(epsilon | cho) is built 
    # in already

    return Kg, Kπprob, v, Q

end


##########################################################################
##########################################################################

function ∑π_ϵ!(muc_ϵ, c, πprob, p, γ)

    Na, Nshocks = size(c)[1:2]

    @inbounds @views for ast = 1:Na

        for shk = 1:Nshocks
            
            muc_ϵ[ast,shk] =  dot((πprob[ast, shk, :]), ( muc.(c[ast, shk, :] , γ) ./ p ))
            # this is the inside part of 54
            # not sure why matmul does not work here...kicks back error

        end
    end

end

##########################################################################
##########################################################################

function make_πprob(vj, σϵ)

    vj_max = maximum(vj, dims = 3) # broadcast, inplace subtract off max value
    # dims = 3 here because countries are stored in 3rd dimension

    foo = vj .- vj_max

    πprob = exp.( foo ./ σϵ ) ./ sum( exp.( foo ./ σϵ ) , dims = 3) 

    return πprob

end

##########################################################################
##########################################################################

function find_asset_position(ga, agrid)

    return findfirst(x -> (x >= ga), agrid)

end

##########################################################################
##########################################################################

function make_Q!(Q, asset_policy, asset_policy_index, πprob, model_params)

    @unpack Na, Nshocks, Ncntry, mc, agrid = model_params
    
    fill!(Q, 0.0) # this is all setup assumeing Q is zero everywehre

    @inbounds @views for cntry = 1:Ncntry

            cntry_counter = Int((cntry - 1)*Nshocks*Na)

            for shk = 1:Nshocks
    
                shk_counter = Int((shk - 1)*Na)

                for ast = 1:Na

                    today = ast + shk_counter + cntry_counter

                    aprime_h = asset_policy_index[ast,shk,cntry]

                    aprime_l = max(aprime_h - 1, 1)
                
                    p = 1.0 - (asset_policy[ast,shk,cntry] - agrid[aprime_l]) / (agrid[aprime_h] - agrid[aprime_l])
                
                    if isnan(p)
                        p = 1.0
                    end

                    for cntryprime = 1:Ncntry

                        cntry_counter_prime = Int((cntryprime - 1)*Nshocks*Na)

                        for shkprime = 1:Nshocks
    
                            shk_counter_prime = Int((shkprime - 1)*Na)

                            for astprime = 1:Na

                                tommorow = astprime + shk_counter_prime + cntry_counter_prime

                                if astprime == aprime_l

                                    Q[today, tommorow] = ( p )*mc.p[shk, shkprime]*πprob[astprime,shkprime,cntryprime]
                                    # given today (a,z,j) = asset choice * prob end up with z' * prob choose varity j'
                                    # or (a',z',j')

                                elseif astprime == aprime_h

                                    Q[today, tommorow] = ( 1.0 - p )*mc.p[shk, shkprime]*πprob[astprime,shkprime,cntryprime]

                                end

                            end 

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

function make_utility!(utility_grid, W, R, p, model_params) 
    # take prices and model parameters and returns utility function
    # R is gross real interest rate ∈ (β, 1/β)
    # W is wage per effeciency unit
    
    @unpack Na, Nshocks, Ncntry, mc, agrid, γ = model_params
    
    a = copy(agrid)
    #assets today
    
    a_prime = transpose(a)
    #assets tomorrow

    shock_level = exp.(mc.state_values)

    for shockstate = 1:Nshocks
    #shock state

        wz = labor_income( shock_level[shockstate] , W)

        for cntry = 1:Ncntry

            c = consumption(R.*a, a_prime, wz, p[cntry])
                # takes assets states, shock state -> consumption from
                # budget constraint
        
            utility_grid[:, :, shockstate, cntry] .= utility.(c, γ)

        end
        
    end
    
end

##########################################################################
##########################################################################

function consumption(Ra, ap, wz, p)
    # constructs consumption via the hh budget constraint
    # .- boadcasting is used to turn vectors (a,ap) into a grid
    # over c (na by na)

    return @. (Ra - ap + wz ) / p
    
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

function bellman_operator_upwind(v, u, mc, β, σϵ) 
    # the value function /bellman operator that takes a v then
    # returns a TV.           
    Na = size(u)[1]
    Nshocks = size(u)[3]
    Ncntry = size(u)[4]
    
    Tv = copy(v)
    # this is here to ensure the output is same type as u...important for
    # autodiff

    Tvj = Array{eltype(u)}(undef, Na, Nshocks, Ncntry)

    for shockstate = 1:Nshocks
        # work through each durable/car state

        βEV = compute_EV(β.*Tv, mc[shockstate, :])
        
        for cntry = 1:Ncntry 
            # work through each shock state

            maximum!(view(Tvj, : , shockstate, cntry), u[:, :, shockstate, cntry] .+ βEV )

        end

        Tvj_max = maximum(Tvj[ :, shockstate, : ], dims = 2)

        foo = Tvj[ :, shockstate, : ] .- Tvj_max

        foo[isnan.(foo)] .= -Inf

        Tv[ :, shockstate] = σϵ.*log.( sum( exp.( ( foo )./ σϵ ) , dims = 2) ) + Tvj_max

    end

    return Tv
      
end
