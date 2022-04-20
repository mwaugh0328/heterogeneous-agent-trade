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

function coleman_operator(policy, R, W, p, model_params)
    # takes gc_j(a,z), v_j(a,z) -> Kgc and Tv
    # multiple dispactch version here for use in fixed point routine

    ######################################################################
    # Organization 

    @unpack agrid, mc, β, γ, σϵ, Na, Nshocks, Ncntry, statesize = model_params

    c = policy[1:Na, :, :]

    v = policy[(Na+1):end, :, :]

    aprime = similar(c)
    Kg = similar(c)
    shocks = exp.(mc.state_values)

    # get choice probabilities
    πprob = make_πprob(v, σϵ)
    muc_ϵ = Array{eltype(R)}(undef, Na, Nshocks)

    ∑π_ϵ!(muc_ϵ, c, πprob, p, γ)
    # this integrates over ϵ

    Emuc = β*R*( matmul( muc_ϵ , mc.p'))
    # this integrates over z

    #Step (3) Work through each county option
    @inbounds @views for cntry = 1:Ncntry

        gc = muc_inverse.( p[cntry] * Emuc, γ)
        # invert from rhs of euler equation
        
        ã = @. (p[cntry] * gc + agrid - W*shocks') / R
        # off budget constraint a = (p_jc_j + a' - w*z ) / R

        # then linear interpolation to get back on grid.
        for shk = 1:Nshocks
    
            foo = LinearInterpolation(ã[:,shk], agrid, extrapolation_bc = (Flat(), Flat()) )

            aprime[:, shk, cntry] = foo.(agrid)
    
            Kg[:, shk, cntry] = @. ( -aprime[:, shk, cntry] + R*agrid + W*shocks[shk] ) / p[cntry]
            # again off budget constraint pc = -a + Ra + wz
    
        end

    end

    # Now I want to infer the value function given updated policy
    Tv = copy(v)

    make_Tv_upwind!(Tv, Kg, aprime, model_params)
    # then Tv = u(g(a,z)) + β*EV

    return vcat(Kg, Tv)

end


##########################################################################
##########################################################################

function coleman_operator(c, v, R, W, p, model_params)
    # Organization 
    @unpack agrid, mc, β, γ, σϵ, Na, Nshocks, Ncntry, statesize = model_params

    aprime = similar(c)
    Kg = similar(c)
    shocks = exp.(mc.state_values)

    # get choice probabilities
    πprob = make_πprob(v, σϵ)
    muc_ϵ = Array{eltype(R)}(undef, Na, Nshocks)

    ∑π_ϵ!(muc_ϵ, c, πprob, p, γ)
    # this integrates over ϵ

    Emuc = β*R*( matmul( muc_ϵ , mc.p'))
    # this integrates over z

    #Step (3) Work through each county option
    @inbounds @views for cntry = 1:Ncntry

        gc = muc_inverse.( p[cntry] * Emuc, γ)
        # invert from rhs of euler equation
        
        ã = @. (p[cntry] * gc + agrid - W*shocks') / R
        # off budget constraint a = (p_jc_j + a' - w*z ) / R

        # then linear interpolation to get back on grid.
        for shk = 1:Nshocks
    
            foo = LinearInterpolation(ã[:,shk], agrid, extrapolation_bc = (Flat(), Flat()) )

            aprime[:, shk, cntry] = foo.(agrid)
    
            Kg[:, shk, cntry] = @. ( -aprime[:, shk, cntry] + R*agrid + W*shocks[shk] ) / p[cntry]
            # again off budget constraint pc = -a + Ra + wz
    
        end

    end

    # Now I want to infer the value function given updated policy
    Tv = copy(v)

    make_Tv_upwind!(Tv, Kg, aprime, model_params)
    # then Tv = u(g(a,z)) + β*EV
    # this function is the bottle neck...worth investing here.
    # why so much memory? 

    return Kg, Tv, aprime

end

##########################################################################
##########################################################################

function make_Tv_upwind!(Tv, Kg, asset_policy, model_params)
    # upwind method that continously updates v as 
    # EV is evaluated....

    @unpack Na, Nshocks, Ncntry, mc, agrid, β, γ, σϵ = model_params

    Ev = 0.0

    @inbounds @views for cntry = 1:Ncntry
        # fix the country

        for shk = 1:Nshocks

            for ast = 1:Na
            
            # here, given aprime, figure out the position
            # appears to by 2x faster not too use function

                aprime_h = searchsortedfirst(agrid, asset_policy[ast, shk, cntry])
            
                aprime_l = max(aprime_h - 1, 1)
            
                p = 1.0 - (asset_policy[ast, shk, cntry] - agrid[aprime_l]) / (agrid[aprime_h] - agrid[aprime_l])
        
                if isnan(p)
                    p = 1.0
                end

            # now work out what happens tomorrow,
            # no need to loop through aprime (we know the position)
            # or country as Ev over variety is the log sum thing

                for shkprime = 1:Nshocks
                    # work through each shock state tommorow to construct
                    # EV
                    
                    Ev += ( p )*mc.p[shk, shkprime]*log_sum_v( Tv[aprime_l, shkprime, :] , σϵ)
                    # so Ev | states today = transition to aprime (p), transition to z',
                    # then multiplies by v tommorow. v tommorow is the log sum thing across different 
                    # options
                
                    Ev += ( 1.0 - p )*mc.p[shk, shkprime]*log_sum_v( Tv[aprime_h, shkprime, :] , σϵ)
                    # note the += here, so we are accumulting this different 
                    # posibilities

                end

                Tv[ast, shk, cntry] = utility_fast(Kg[ast, shk, cntry], γ) + β*Ev

                #Then the vj = uj + βEV

                Ev = 0.0

            end

        end

    end

end

##########################################################################
##########################################################################



function make_Q!(Q, household, model_params)

    @unpack asset_policy, πprob = household
    @unpack Na, Nshocks, Ncntry, mc, agrid = model_params

    fill!(Q, 0.0) # this is all setup assumeing Q is zero everywehre
    # Q is size Na X Nshocks. Country variety not being tracked.

    for cntry = 1:Ncntry
        # fix a country and then work through each shock, asset situation

        for shk = 1:Nshocks

            shk_counter = Int((shk - 1)*Na)

            for ast = 1:Na

                today = ast + shk_counter

                aprime_h = searchsortedfirst(agrid, asset_policy[ast, shk, cntry])
                    # asset choice for ast, shk, and variety choice
            
                aprime_l = max(aprime_h - 1, 1)
                
                p = 1.0 - (asset_policy[ast, shk, cntry] - agrid[aprime_l]) / (agrid[aprime_h] - agrid[aprime_l])

                if isnan(p)
                    p = 1.0
                end

                for shkprime = 1:Nshocks
                        # then work through tomorrow 

                    shk_counter_prime = Int((shkprime - 1)*Na)

                    tommorow = aprime_l + shk_counter_prime

                    Q[today, tommorow] += ( p )*mc.p[shk, shkprime]*πprob[ast, shk, cntry]
                        # given today (a,z,j) = asset choice(a,z,j) * prob end up with z' * prob choose varity j
                        # given today (a,z), then the += accumulation here picks up as we work through cntry

                    tommorow = aprime_h + shk_counter_prime

                    Q[today, tommorow] += ( 1.0 - p )*mc.p[shk, shkprime]*πprob[ast, shk, cntry]


                end

            end

        end
    end

end

##########################################################################
##########################################################################


function make_πprob(vj, σϵ)
    
        foo = vj .- maximum(vj, dims = 3)

        return exp.( foo ./ σϵ ) ./ sum( exp.( foo ./ σϵ ) , dims = 3) 
    
end

##########################################################################
##########################################################################

function log_sum_v(vj, σϵ)

    vj_max = maximum(vj)

    foo = vj .- vj_max

    return σϵ*log( sum( exp.( ( foo )/ σϵ )) ) + vj_max

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

function utility_fast(c, γ)
    # maps consumption into utility with the CRRA specification
    # log it γ is close to one

    if γ ≈ 1.0
        
        @fastmath log(c) 

    else
        @fastmath c^( one(γ) - γ) / (one(γ) - γ)

    end

end


##########################################################################
##########################################################################

function utility(c, γ)
    # maps consumption into utility with the CRRA specification
    # log it γ is close to one

    if γ ≈ 1.0
        
        (c < 1e-10 ? -Inf : log(c) )

    else
        (c < 1e-10 ? -Inf : c^( 1.0 - γ) / (1.0 - γ))

    end

end

##########################################################################
##########################################################################

function make_utility!(utility_grid, R, W, p, model_params) 
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

function log_sum_column(vj, σϵ)

    Na = size(vj)[1]
    Nshocks = size(vj)[2]

    vj_max = maximum(vj, dims = 3)

    foo = vj .- vj_max

    return reshape(σϵ*log.( sum( exp.( ( foo ) / σϵ ) , dims = 3) ) + vj_max, Na, Nshocks)

end

##########################################################################
##########################################################################

function make_state_index!(state_index, model_params)
    
    @unpack Na, Nshocks = model_params

    for shk = 1:Nshocks

        shk_counter = Int((shk - 1)*Na)

        for ast = 1:Na

            today = ast + shk_counter

            state_index[today] = (ast,shk)

        end

    end

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
