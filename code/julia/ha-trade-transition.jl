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

function transition_path(xxx, Rpath, τrsf_path, trp_values, foo_params; display = false)
    # multiple dispatch version for use in solver

    #####################################################################################################
    # ORGANIZATION NEED TO BE FIXED

    @unpack Ngoods, Na, Nshocks, ς = foo_params
    @unpack hh_end, dist₀, Rend, W, T, τ, τrsf, G, pend = trp_values

    # R = vcat([R₀], Rpath, [Rend])
    # πrft = xxx[1:(T)]
    # p = reshape(xxx[T + 1 : end], Ngoods, T)
    # this is for situation with initial pinned down

    R = vcat(Rpath, [Rend])# add the final period
    πrft = xxx[1:T]
    p = reshape(xxx[T + 1 : end], Ngoods, T)
    p = hcat(p, pend) # add the final period

    τ = τ.*ones(T+1)
    τrsf = τrsf_path
    # this is status quo transfers + new transfers on the path 

    @assert length(R) ≈ T + 1
    @assert length(πrft) ≈ T
    @assert size(p)[2] ≈ T + 1 #need to have end point here

    hh = Array{household{eltype(W)}}(undef, T+1, Ncntry) #### LOOK AT THIS
    # add in the end period as the + 1

    p̂ = Array{eltype(W)}(undef, Ngoods, T)

    net_asset_demand = Array{eltype(W)}(undef, T)

    p_diff = Array{eltype(W)}(undef, Ngoods, T)

    πrft_diff = Array{eltype(W)}(undef, T)

    B = Array{eltype(W)}(undef, T+1)
    B[1] = foo_params.B

    for cntry = 1:Ncntry # for each country

        hh[end, cntry] = hh_end[cntry]
    # this is the household at the end

        λ[1, cntry] = deepcopy(dist₀[cntry].λ)

    end

    Q = Array{Array{Float64,2}}(undef, T, Ncntry)
    #Q = Array{Float64}(undef, Na*Nshocks, Na*Nshocks)

  for fwdate = 1:T

    #### Country dimension needs to be changed

        Q[fwdate, cntry] = Array{Float64}(undef, Na*Nshocks, Na*Nshocks, Ncntry)

    end


    #####################################################################################################
    # This is the backward step: sovle hh problem at T then use colman operator to work backwards

    for bwdate = (T):-1:1 # do this for each T

        for cntry = 1:Ncntry # for each country

            p[:, bwdate] = make_p(W[:, bwdate], TFP[:, bwdate], d[cntry, :, bwdate], tariff[cntry, :, bwdate] ) # need T add t dimension
    
            ψ = make_ψ(cntry, ψslope.*TFP[cntry].^(1.0 - γ), hh_params) # need T add t dimension
            # this creates the z quality shifter
            # scaled in a way that is invariant to level of TFP
    
            agrid = make_agrid(hh_params, TFP[cntry]) # need T add t dimension
            # this creates teh asset grid so it's alwasy a fraction of home labor income
    
            foo_hh_params = household_params(hh_params, agrid = agrid, 
                    TFP = TFP[cntry], L = L[cntry], σϵ = σϵ*(TFP[cntry]^(1.0 - γ)), ψ = ψ) # need T add t dimension

            #### THIS IS WHERE OUR NEW ONE STEP WOULD GO
            hh[bwdate] = one_step_itteration(hh[bwdate + 1].consumption_policy, hh[bwdate + 1].Tv, # consumption, values at date t+1
                    R[bwdate], R[bwdate + 1], # returns at date t and t + 1
                    W, πrft[bwdate], # factor prices at date t
                    p[:, bwdate] , p[:, bwdate + 1], foo_hh_params) # goods prices at date t and t+1

            ### THIS WOULD NEED TO HAVE COUNTRY DIMENSION
    
        end
        
    end

    # this constructs the transition matrix, here there is no 
    # recursive relationship, so it can be multi-threaded
    for fwdate = 1:T

        foo_model_params = model_params(foo_params, τrsf = τrsf[fwdate], τ = τ[fwdate])

        make_Q!(Q[fwdate], hh[fwdate], foo_model_params) # THIS WOULD NEED TO BE BY COUNTRY

    end

    #####################################################################################################
    # This is the forward step, so given an initial distribution, take hh decicion rules and push forward

    for fwdate = 1:T
        # so when date > T as we run it out, just grab stuff from end in policy functions or parameter

        for cntry = 1:Ncntry

            p = make_p(W, TFP, d[cntry, :], tariff[cntry, :] ) #ADD T dimension
    
            ψ = make_ψ(cntry, ψslope.*TFP[cntry].^(1.0 - γ), hh_params) #ADD T dimension
    
            agrid = make_agrid(hh_params, TFP[cntry]) #ADD T dimension
    
            foo_hh_params = household_params(hh_params, agrid = agrid, 
                    TFP = TFP[cntry], L = L[cntry], σϵ = σϵ*(TFP[cntry]^(1.0 - γ)), ψ = ψ) #ADD T dimension
    
            output, tradestats = aggregate(R[cntry], W[cntry], p, τ[cntry], tariff, cntry, 
                hh[cntry, fwdate], distribution(Q[cntry, fwdate], λ[cntry, fwdate], dist₀.state_index), foo_hh_params)
            #ADD T dimension
    
            Y[cntry, fwdate] = output.production
    
            tradeflows[cntry, :, fwdate] = tradestats.bilateral_imports
        
            A_demand[cntry, fwdate] = output.Aprime

            λ[cntry, fwdate + 1] .= law_of_motion(λ[cntry, fwdate] , transpose(Q[fwdate]))
    
        end            
            #then push forward

        goods_market[:, fwdate] = Y[:, fwdate] .- vec(sum(tradeflows[:, :, fwdate] , dims = 1))
        
        asset_market[fwdate] = sum(A_demand[:, fwdate])

    end

    return goods_market,  asset_market # the goods market should by Ncntry by T, asset market should be by T 

end

#####################################################################################################