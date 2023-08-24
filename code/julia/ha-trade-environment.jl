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
#Pkg.pin(name = "Interpolations", version = "0.13.1")
using LinearAlgebra
using LoopVectorization
using FiniteDifferences
using Random
using DelimitedFiles

include("mix-MarkovChain.jl")

##########################################################################
@with_kw struct household_params
    TFP::Float64 = 1.0
    L::Float64 = 1.0
    β::Float64 = 0.95
    γ::Float64 = 2.0
    ϕ::Float64 = 0.0
    amax::Float64 = 8.0
    Ncntry::Int64 = 2
    σϵ::Float64 = 0.25
    Na::Int64 = 50
    agrid::Array{Float64, 1} = convert(Array{Float64, 1}, range(-ϕ*TFP, amax*TFP, length = Na))
    Nar::Int64 = 5
    Nma::Int64 = 2
    Nshocks::Int64 = Nar*Nma
    statesize::Int64 = Int(Na*Nshocks*Ncntry)
    ρ::Float64 = 0.90
    σar::Float64 = 0.039^(0.5)
    σma::Float64 = 0.0522^(0.5)
    mc::MarkovChain{Float64, Matrix{Float64}, Vector{Float64}} = mMarkovChain(Nar,Nma,ρ,σar,σma)
    ψ::Array{Float64, 3} = zeros(Na, Nshocks, Ncntry)
    ψslope::Float64 = 0.0
    λτ::Float64 = 1.0
end

@with_kw struct country_params
    Ncntry::Int64 = 2
    TFP::Array{Float64, 1} = ones(Ncntry)
    L::Array{Float64, 1} = ones(Ncntry)
    d::Array{Float64, 2} = ones(Ncntry,Ncntry)
    tariff::Array{Float64, 2} = zeros(Ncntry,Ncntry)
 end

##########################################################################

function coleman_operator(policy, R, W, p, τ, model_params)
    # multiple dispatch version that directly takes policy functions
    # single R so this means assumption is a stationary setting

    c = policy[1:model_params.Na, :, :]

    v = policy[(model_params.Na+1):end, :, :]

    Kg, Tv = coleman_operator(c, v, R, W, p, τ, model_params)[1:2]

    return vcat(Kg, Tv)

end

##########################################################################
##########################################################################

function coleman_operator(c, v, R, W, p, τ, model_params)
    # Organization 
    @unpack agrid, mc, β, γ, σϵ, ψ, λτ, Na, Nshocks, Ncntry = model_params

    aprime = similar(c)
    Kg = similar(c)
    shocks = exp.(mc.state_values)

    muc_ϵ = Array{eltype(R)}(undef, Na, Nshocks)

    ã = Array{eltype(R)}(undef, Na, Nshocks)
    gc = Array{eltype(c)}(undef, Na, Nshocks)

    Emuc = Array{eltype(R)}(undef, Na, Nshocks)

    # get choice probabilities

    # if sum(isnan.(v)) != 0
    #     println(v)
    # end

    πprob = make_πprob(v, σϵ, ψ)

    ∑π_ϵ!(muc_ϵ, c, πprob, p, γ)
    # this integrates over ϵ

    Emuc .= β*R*λτ*( matmul( muc_ϵ , mc.p') )
    # λτ is porportional tax / subsidy on 
    # all income. used for welfare analysis

    #Step (3) Work through each county option
     @inbounds @views for cntry = 1:Ncntry

        muc_inverse!( gc, p[cntry] * Emuc, γ)

        ã .= @. (p[cntry] * gc + agrid - λτ*W*shocks' - τ) / ( R*λτ )
        # off budget constraint a = (p_jc_j + a' - w*z - τ ) / R

        # then linear interpolation to get back on grid.
        for shk = 1:Nshocks

            # if issorted(ã[:, shk, cntry]) == false
            #     println(c[1,1,1])
            #     println(πprob[1,1,1])
            #     println(v[1,1,1])
            # end

            foo = LinearInterpolation(ã[:, shk], agrid, extrapolation_bc = (Flat(), Flat()) )
    
            aprime[:, shk, cntry ] .= foo.(agrid)

            Kg[:, shk, cntry] .= @. ( -aprime[:, shk, cntry] + λτ*(R*agrid + W*shocks[shk]) + τ ) / p[cntry]
            # again off budget constraint pc = -a' + Ra + wz + τ 
    
        end

    end

    # Now I want to infer the value function given updated policy
    Tv = similar(v)

    make_Tv!(Tv, v, Kg, aprime, πprob, 1.0, ψ, model_params)
    # then Tv = u(g(a,z)) + β*EV
    # this function is the bottle neck...worth investing here.
    # why so much memory? 

    return Kg, Tv, aprime

end

##########################################################################
##########################################################################

function make_Tv!(Tv, v, Kg, asset_policy, πprob, model_params)
    # multiple dispatch for no quality case

    make_Tv!(Tv, v, Kg, asset_policy, πprob, 1.0, 0.0, model_params)

end

function make_Tv!(Tv, v, Kg, asset_policy, πprob, λ, ψ, model_params)
    # constructs the choice specific value functions

    @unpack Na, Nshocks, Ncntry, mc, agrid, β, γ, σϵ = model_params

    Ev = Array{eltype(v)}(undef, Ncntry) 
    # need to have it like this to mulithread
    fill!(Ev, 0.0)

    Φ = reshape(alt_log_sum(πprob, v, σϵ, ψ), Na, Nshocks)

    @inbounds @views for cntry = 1:Ncntry
        # fix the country

        for shk = 1:Nshocks

            for ast = 1:Na
            
            # here, given aprime, figure out the position
            # appears to by 2x faster not too use function

                aprime_h = searchsortedfirst(agrid, asset_policy[ast, shk, cntry])
                #searchsortedfirst.(Ref(agrid), asset_policy[:,shk, cntry])
                # broadcaseted version
            
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
                    prob_eprime = mc.p[shk, shkprime]

                    Ev[cntry] += ( p )*prob_eprime*Φ[aprime_l, shkprime]
                    # so Ev | states today = transition to aprime (p), transition to z',
                    # then multiplies by v tommorow. v tommorow is the log sum thing across different 
                    # options
                
                    Ev[cntry] += ( 1.0 - p )*prob_eprime*Φ[aprime_h, shkprime]
                    # note the += here, so we are accumulting this different 
                    # posibilities

                end

                Tv[ast, shk, cntry] = utility(λ*Kg[ast, shk, cntry], γ) + β*Ev[cntry]

                #Then the vj = uj + βEV
                #This does not include quality ψ

                Ev[cntry] = 0.0

            end

        end

    end

end

##########################################################################
##########################################################################
function log_sum_v(vj, σϵ, Ncntry)
    # this just does a loop on it.
    # faster and way more memory effecient

    foo = 0.0

    vj_max = maximum(vj)
    # this part slows down by 2X...
    # is it necessary? one thought is make
    # it relative to home country

    @inbounds @turbo for xxx = 1:Ncntry

        foo += exp( ( vj[xxx] - vj_max ) / σϵ )

    end

    return σϵ*log( foo ) + vj_max

end

function log_sum_v(vj, σϵ)

    #foo = @. exp( ( vj ) / σϵ )
    # not haveing ./ on the σϵ was a problem

    vj_max = maximum(vj, dims = 3)

    foo = vj .- vj_max

    return σϵ*log.( sum( exp.( ( foo ) / σϵ ) , dims = 3) ) .+ vj_max
    # this is the prolem code...use alt_log_sum below
end

function alt_log_sum(household, σϵ)

    @unpack πprob, Tv = household

    foo = πprob .* ( Tv .- σϵ.*log.(πprob))

    return sum(foo, dims = 3)

end

function alt_log_sum(πprob, Tv, σϵ)

    return alt_log_sum(πprob, Tv, σϵ, 0.0)

end

function alt_log_sum(πprob, Tv, σϵ, ψ)

    foo = πprob .* ( ( ψ .+ Tv ) .- σϵ.*log.(πprob))

    return sum(foo, dims = 3)

end


##########################################################################
##########################################################################



function make_Q!(Q, household, model_params)

    @unpack asset_policy, πprob = household
    @unpack Na, Nshocks, Ncntry, mc, agrid = model_params

    fill!(Q, 0.0) # this is all setup assumeing Q is zero everywehre
    # Q is size Na X Nshocks. Country variety not being tracked.

    @inbounds @views for cntry = 1:Ncntry
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

                π_a_z_j = πprob[ast, shk, cntry]

                for shkprime = 1:Nshocks
                        # then work through tomorrow
                        
                    prob_zprime_j = mc.p[shk, shkprime]*π_a_z_j   

                    shk_counter_prime = Int((shkprime - 1)*Na)

                    tommorow = aprime_l + shk_counter_prime

                    Q[today, tommorow] += ( p )*prob_zprime_j 
                        # given today (a,z,j) = asset choice(a,z,j) * prob end up with z' * prob choose varity j
                        # given today (a,z), then the += accumulation here picks up as we work through cntry

                    tommorow = aprime_h + shk_counter_prime

                    Q[today, tommorow] += ( 1.0 - p )*prob_zprime_j 


                end

            end

        end
    end

end

##########################################################################
##########################################################################

# function make_πprob(vj, σϵ)

#     # if sum(isnan.(vj)) != 0
#     #     println("nan showing up in make prob")
#     # end

#     # .- maximum(vj, dims = 3)
#     # this is strange, was causing a nan to show up
#     # can't replicate this behavior in simple examples.

#     foo = vj .- maximum(vj, dims = 3)

#     foo .= exp.( foo / σϵ)

#     # if sum(isnan.(foo)) != 0
#     #     println(vj[1,1,:])
#     #     println(foo[1,1,:])
#     #     println("nan showing up")
#     # end

#     # goof = sum( foo, dims = 3)

#     # foobar = foo ./ goof

#     # if sum(isnan.(foobar)) != 0
#     #     println(vj[1,1,:])
#     #     println(foo[1,1,:])
#     #     println("nan in denomenator showing up")
#     # end

#     return foo ./ sum( foo, dims = 3)
   
# end

function make_πprob(vj, σϵ)

    return make_πprob(vj, σϵ, 0.0)

end

function make_πprob(vj, σϵ, ψ)
    # constructs πs with quality shifter
    # exp ( (ψ + v(a,z,j) ) / σϵ) / Φ(a,z)

    foo = (ψ .+ vj) .- maximum( (ψ .+ vj) , dims = 3)

    foo .= exp.( foo / σϵ)

    return foo ./ sum( foo, dims = 3)
   
end

##########################################################################
##########################################################################

function make_ψ(home, ψslope, model_params)

    @unpack Na, Nar, Nma, ρ, Nshocks, Ncntry, σar, = model_params

    ~ , perm_z = rouwenhorst(Nar, ρ, σar, 0.0)
    
    ψ = similar(model_params.ψ)
    fill!(ψ, 0.0)

    perm_z = repeat(perm_z, inner = Nma, outer = 1)
    
    slope = ψslope.*perm_z

    @assert length(slope) ≈ Nshocks
    
    for cntry = 1:Ncntry

        if home == cntry
        
            for ast = 1:Na

                for shk = 1:Nshocks

                    ψ[ast, shk, cntry] = -0.0 + slope[shk]

                end
            end
        end
    end

    return ψ

end

##########################################################################
##########################################################################

function ∑π_ϵ!(muc_ϵ, c, πprob, p, γ)

    Na, Nshocks = size(c)[1:2]

    @inbounds @views for ast = 1:Na

        for shk = 1:Nshocks
            
            #muc_ϵ[ast,shk] = dot((πprob[ast, shk, :]), ( muc.(c[ast, shk, :] , γ) ./ p ))
            muc_ϵ[ast,shk] = mydot(πprob[ast, shk, :], c[ast, shk, :], p, γ) 

            # this is the inside part of 54
            # not sure why matmul does not work here...kicks back error

        end
    end

end

function mydot(πprob, c, p, γ)
    s = 0.0

    @turbo for i ∈ eachindex(πprob)

        s += πprob[i] * ( muc(c[i], γ) / p[i] )

    end
    s
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

    @fastmath @. c^(-γ) 

end

##########################################################################
##########################################################################
function muc_inverse!(out, c, γ)
    # the marginal utility of consumption with
   # CRRA preferences

   out .= @. c^( -one(γ) / γ )

end

function muc_inverse(c, γ)
     # the marginal utility of consumption with
    # CRRA preferences

    @. c^( 1.0 / -γ )

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
        
        (c ≤ 0.0 ? -Inf : log(c) )

    else
        (c ≤ 0.0 ? -Inf : c^( one(γ) - γ)  / (one(γ) - γ))

    end

end

##########################################################################
##########################################################################

function make_utility!(utility_grid, R, W, p, τ,  model_params) 
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

    return @. (Ra - ap + wz) / p
    
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

    #foo = vj .- vj_max

    return reshape(σϵ*log.( sum( exp.( ( vj .- vj_max ) / σϵ ) , dims = 3) ) + vj_max, Na, Nshocks)

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

##########################################################################
##########################################################################

function make_d!(d, dij)
    Ncntry = size(d)[1]

    for cntry = 1:Ncntry

        d[1:end .!= cntry, :] .= dij

    end

end

function make_p(W, TFP, d, tariff)

    return (W ./ TFP) .* d .* (1.0 .+ tariff)

end


##########################################################################

function mean_z(model_params)
    # gets the mean value from the z-shocks
    # from the amador-allerano code

    @unpack Nshocks, mc = model_params

    shocks = exp.(mc.state_values)

    invZ = real(inv(eigvecs(mc.p))[Nshocks,:])
    invZ .= invZ ./ sum(invZ)

    return invZ' * shocks

end