struct trans_path_values
    hh_end::Array{household{Float64}, 1} # check this later
    dist₀::Array{distribution{Float64}, 1} # check this later
    R₀::Array{Float64, 1} # Ncntry by 1
    Rend::Array{Float64, 1} # Ncntry by 1
    Wend::Array{Float64, 1} # Ncntry by 1
    T::Int64
    τ::Array{Float64, 1} # this is transfer, Ncntry by 1, set to [0.0; 0.0] for now
end

#####################################################################################################

function push_foward!(λ, Q, household, model_params)
    # Pushes the economy forward
    # (1) Take policy functions and a distribution L -> aggregates today.
    # Policy functions -> Transition probability Q
    # Distribution today + Q -> Distribution tomorrow. 
    make_Q!(Q, household, model_params)

    # this had an alt_make_Q, not sure what difference is, more complicated
    # 
    
    # Then push the distribution forward
    λ = law_of_motion(λ , transpose(Q))
    
end

#####################################################################################################

function one_step_itteration(cₜ₊₁, vₜ₊₁, Rₜ, Rₜ₊₁, Wₜ, pₜ, pₜ₊₁, τ, model_params)
    # used to work backward. Give me a policy function and V at date
    # t + 1, I return a policy function and V for t
    
    Kgcₜ , Tvₜ , Kgaₜ = coleman_operator(cₜ₊₁, vₜ₊₁, Rₜ, Rₜ₊₁, Wₜ, pₜ, pₜ₊₁, τ, model_params)

    πprob = make_πprob(Tvₜ , model_params.σϵ, model_params.ψ)

    return household(Kgaₜ , Kgcₜ , πprob, Tvₜ)

end

#####################################################################################################
# here the idea is there is a change in trade costs, that is the exogenous path changed
# it would be set up so d path, W path, and R path that goes, so set up argumetns that fixes that

function transition_path(xxx, d_path, trp_values, hh_params, cntry_params; display = false)
    # multiple dispatch version to deal with endogenous path of R 

    @unpack T = trp_values

    Rpath = xxx[1:(T-1)]

    goods_market, asset_market = transition_path(xxx[T:end], Rpath, d_path, trp_values, hh_params, cntry_params; display = display)
    
    return vcat(goods_market, asset_market[1:(T - 1)])

end

#####################################################################################################

function transition_path_only_assetmarket(xxx, w_path, d_path, trp_values, hh_params, cntry_params; display = false)
    # multiple dispatch version to deal with endogenous path of R 
    # Rpath = xxx[1:T]

    asset_market = transition_path(w_path, xxx, d_path, trp_values, hh_params, cntry_params; display = display)[2]
    
    return asset_market[1:(trp_values.T - 1)]

end

#####################################################################################################
#####################################################################################################

function transition_path(xxx, Rpath, d_path, trp_values, hh_params, cntry_params; display = false)
    # multiple dispatch version for use in solver

    @unpack ψslope, γ, σϵ, Ncntry, Na, Nshocks, mc = hh_params
    @unpack hh_end, dist₀, R₀, Rend, Wend, T, τ = trp_values
    @unpack TFP, L, tariff = cntry_params # decide if want to put TFP and L in 'trp_values'

    # R = vcat([R₀], Rpath, [Rend])
    # πrft = xxx[1:(T)]
    # p = reshape(xxx[T + 1 : end], Ngoods, T)
    # this is for situation with initial pinned down

    R = reshape(Rpath, Ncntry, T-1)
    R = hcat(R₀, Rpath, Rend) # add the final period
    W = reshape(xxx[:], Ncntry, T) # we are finding W path such that markets clear at all date, here we feed in xxx as 2T by 1 vector
    W = hcat(W, Wend)

    @assert length(R[1,:]) ≈ T + 1
    @assert length(W[1,:]) ≈ T + 1 

    TFP = TFP.*ones(Ncntry, T+1)
    τ = τ.*ones(Ncntry, T+1) # transfer
    L = L.*ones(Ncntry, T+1)
    tariff = tariff.*ones(Ncntry, Ncntry, T+1)

    hh = Array{household{eltype(W)}}(undef, Ncntry, T+1) #### LOOK AT THIS
    # add in the end period as the + 1

    λ = Array{eltype(W)}(undef, Na*Nshocks, Ncntry, T+1) # This is different from PI, here we define λ as a 3-dimensional object

    # Aggregate stuff holders
    Y = Array{eltype(W)}(undef, Ncntry, T)
    tradeflows = Array{eltype(W)}(undef, Ncntry, Ncntry, T)
    A_demand = Array{eltype(W)}(undef, Ncntry, T)
    goods_market = Array{eltype(W)}(undef, Ncntry, T)
    asset_market = Array{eltype(W)}(undef, T)

    @views @inbounds for cntry = 1:Ncntry # for each country

        hh[cntry, end] = hh_end[cntry]
        # this is the household at the end

        λ[:, cntry, 1] .= deepcopy(dist₀[cntry].λ)
        # this is dist. at beginning

    end

    #####################################################################################################
    # This is the backward step: solve hh problem at T then use colman operator to work backwards

    @time @views @inbounds for bwdate = (T):-1:1 # do this for each T

        for cntry = 1:Ncntry # for each country

            pₜ = make_p(W[:, bwdate], TFP[:, bwdate], d_path[cntry, :, bwdate], tariff[cntry, :, bwdate]) # need T add t dimension
    
            pₜ₊₁ = make_p(W[:, bwdate + 1], TFP[:, bwdate + 1], d_path[cntry, :, bwdate + 1], tariff[cntry, :, bwdate + 1])
    
            foo_hh_params = household_params(hh_params, agrid = make_agrid(hh_params, TFP[cntry, bwdate]), 
                    TFP = TFP[cntry, bwdate], L = L[cntry, bwdate], σϵ = σϵ*(TFP[cntry, bwdate]^(1.0 - γ)),
                     ψ = make_ψ(cntry, ψslope.*TFP[cntry, bwdate].^(1.0 - γ), hh_params) ) # need T add t dimension

            hh[cntry, bwdate] = one_step_itteration(hh[cntry, bwdate + 1].cons_policy, hh[cntry, bwdate + 1].Tv, # consumption, values at date t+1
                    R[cntry, bwdate], R[cntry, bwdate + 1], # returns at date t and t + 1
                    W[cntry, bwdate], # factor prices at date t
                    pₜ , pₜ₊₁, τ[cntry, bwdate], foo_hh_params) # goods prices at date t and t+1
        end
        
    end

    # this constructs the transition matrix, here there is no 
    # recursive relationship, so it can be multi-threaded

    Q = Array{Float64}(undef, Na*Nshocks, Na*Nshocks)

    @time @views @inbounds for fwdate = 1:T

        for cntry = 1:Ncntry

            foo_hh_params = household_params(hh_params, agrid = make_agrid(hh_params, TFP[cntry, fwdate]), 
                        TFP = TFP[cntry, fwdate], L = L[cntry, fwdate], σϵ = σϵ*(TFP[cntry, fwdate]^(1.0 - γ)),
                         ψ = make_ψ(cntry, ψslope.*TFP[cntry, fwdate].^(1.0 - γ), hh_params))

            make_Q!(Q, hh[cntry, fwdate], foo_hh_params) 

            λ[:, cntry, fwdate + 1] .= law_of_motion(λ[:, cntry, fwdate] , transpose(Q) )
            
        end

    end

    #####################################################################################################
    # This is the forward step, so given an initial distribution, take hh decision rules and push forward

    @time @inbounds @views for fwdate = 1:T
        # so when date > T as we run it out, just grab stuff from end in policy functions or parameter

        for cntry = 1:Ncntry

            pₜ = make_p(W[:, fwdate], TFP[:, fwdate], d_path[cntry, :, fwdate], tariff[cntry, :, fwdate])
    
            foo_hh_params = household_params(hh_params, agrid = make_agrid(hh_params, TFP[cntry, fwdate]), 
                        TFP = TFP[cntry, fwdate], L = L[cntry, fwdate], σϵ = σϵ*(TFP[cntry, fwdate]^(1.0 - γ)),
                         ψ = make_ψ(cntry, ψslope.*TFP[cntry, fwdate].^(1.0 - γ), hh_params))
    
            output, tradestats = aggregate(R[cntry, fwdate], W[cntry, fwdate], pₜ, τ[cntry, fwdate], tariff[:,:,fwdate], cntry, 
                hh[cntry, fwdate], distribution(Q, λ[:, cntry, fwdate], dist₀[cntry].state_index), foo_hh_params)

            # this is the next step to simplify...

            Y[cntry, fwdate] = output.production
    
            tradeflows[cntry, :, fwdate] = tradestats.bilateral_imports
        
            A_demand[cntry, fwdate] = output.Aprime
    
        end            
            #then push forward

        goods_market[:, fwdate] .= Y[:, fwdate] .- vec(sum(tradeflows[:, :, fwdate] , dims = 1))
        
        asset_market[fwdate] = sum(A_demand[:, fwdate])

    end

    return goods_market,  asset_market # the goods market should by Ncntry by T, asset market should be by T 

end

#####################################################################################################