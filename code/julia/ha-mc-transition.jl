struct trans_path_values
    hh_end::household{Float64}  
    dist₀::distribution{Float64} 
    R₀::Float64 
    Rend::Float64
    W::Float64 
    T::Int64
    τ::Float64
    τrsf::Float64
    G::Float64 
    pend::Array{Float64} 
end

struct hhmuc
    poor::Float64
    rich::Float64
    middle::Float64
end


#####################################################################################################

function transition_path(xxx, τrsf_path, trp_values, foo_params; display = false)
    # multiple dispatch version to deal with endogenous path of R 
    Rpath = xxx[1:T]

    out1, out2 = transition_path(xxx[T+1:end], Rpath, τrsf_path, trp_values, foo_params; display = display)
    
    return vcat(out1, out2)

end

#####################################################################################################

function transition_path(xxx, Rpath, τrsf_path, trp_values, foo_params; display = false)
    # multiple dispatch version for use in solver

    #####################################################################################################
    # This is just getting organized

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

    hh = Array{household{eltype(W)}}(undef, T+1, 1)
    # add in the end period as the + 1

    p̂ = Array{eltype(W)}(undef, Ngoods, T)

    net_asset_demand = Array{eltype(W)}(undef, T)

    p_diff = Array{eltype(W)}(undef, Ngoods, T)

    πrft_diff = Array{eltype(W)}(undef, T)

    B = Array{eltype(W)}(undef, T+1)
    B[1] = foo_params.B

    hh[end] = hh_end
    # this is the household at the end

    λ = deepcopy(dist₀.λ)
    Q = Array{Array{Float64,2}}(undef, T)
    #Q = Array{Float64}(undef, Na*Nshocks, Na*Nshocks)

    @inbounds for fwdate = 1:T

        Q[fwdate] = Array{Float64}(undef, Na*Nshocks, Na*Nshocks)

        B[fwdate + 1] = ( τ[fwdate].*W - τrsf[fwdate] ) - G + R[fwdate]*(B[fwdate])
                
        τ[fwdate + 1] = τ[1]*( B[fwdate+1] / B[1] )^ς

    end

    # println(B)
    # println(τ)

    Aprime = Array{Float64}(undef, T, 1)
    πrft_hat = similar(Aprime)

    #####################################################################################################
    # This is the backward step: sovle hh problem at T then use colman operator to work backwards

    @inbounds for bwdate = (T):-1:1

        foo_model_params = model_params(foo_params, τrsf = τrsf[bwdate], τ = τ[bwdate]) # transfer showing up at date t

        # one_step_itteration(cₜ₊₁, vₜ₊₁, Rₜ, Rₜ₊₁, W, πrft, pₜ, pₜ₊₁, model_params) # this is the funciton
        
        hh[bwdate] = one_step_itteration(hh[bwdate + 1].consumption_policy, hh[bwdate + 1].Tv, # consumption, values at date t+1
                    R[bwdate], R[bwdate + 1], # returns at date t and t + 1
                    W, πrft[bwdate], # factor prices at date t
                    p[:, bwdate] , p[:, bwdate + 1], foo_model_params) # goods prices at date t and t+1

        # then this works backward, so it delivers household policy functions at date t
        
    end

    # this constructs the transition matrix, here there is no 
    # recursive relationship, so it can be multi-threaded

    @inbounds Threads.@threads for fwdate = 1:T

        foo_model_params = model_params(foo_params, τrsf = τrsf[fwdate], τ = τ[fwdate])

        make_Q!(Q[fwdate], hh[fwdate], foo_model_params)

    end

    #####################################################################################################
    # This is the forward step, so given an initial distribution, take hh decicion rules and push forward

    if display == false
        
        @inbounds for fwdate = 1:T
        # so when date > T as we run it out, just grab stuff from end in policy functions or parameter

            foo_model_params = model_params(foo_params, τrsf = τrsf[fwdate], τ = τ[fwdate])

            p̂[:, fwdate], Aprime[fwdate], πrft_hat[fwdate] = aggregate_lite(W, p[:, fwdate], hh[fwdate], 
            distribution(Q[fwdate], λ, dist₀.state_index), foo_model_params)

            #somehow this is a bottle neck
        
            net_asset_demand[fwdate] = Aprime[fwdate] + B[fwdate + 1]
            # 

            @views p_diff[:, fwdate] = p[:, fwdate] .- p̂[:, fwdate]
                
            πrft_diff[fwdate] = πrft[fwdate] - πrft_hat[fwdate]

            λ .= law_of_motion(λ , transpose(Q[fwdate]))
            #then push forward

    end

        return vcat(πrft_diff, p_diff[:]),  net_asset_demand[1:T]

elseif display == true

        output = Array{NIPA}(undef, T, 1)
        hh_muc = Array{hhmuc}(undef, T, 1)

        @inbounds for fwdate = 1:T
            # so when date > T as we run it out, just grab stuff from end in policy functions or parameter
            
            #alt_make_Q!(Q, hh[fwdate], model_params)
    
            foo_model_params = model_params(foo_params, τrsf = τrsf[fwdate], τ = τ[fwdate])
    
    
            output[fwdate] = aggregate(R[fwdate], W, πrft[fwdate], p[:, fwdate], hh[fwdate], 
                        distribution(Q[fwdate], λ, dist₀.state_index), foo_model_params)[2]

            df_firm = make_firm_dataframe(output[fwdate], p[:, fwdate], foo_model_params)

            df = alt_hh_dataframe(hh[fwdate], distribution(Q[fwdate], λ, dist₀.state_index), 
                df_firm, R[fwdate], W, πrft[fwdate], foo_model_params)

            poor, rich, middle = muc_stats(df)

            hh_muc[fwdate] = hhmuc(poor, rich, middle)

            λ .= law_of_motion(λ , transpose(Q[fwdate]))
        end

        return output, B, hh_muc

    end

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

function one_step_itteration(cₜ₊₁, vₜ₊₁, Rₜ, Rₜ₊₁, W, πrft, pₜ, pₜ₊₁, model_params)
    # used to work backward. Give me a policy vunction and V at date
    # t + 1, I return a policy function and V for t
    
    Kgcₜ , Tvₜ , Kgaₜ , Emucₜ₊₁ , V̄ₜ₊₁ = coleman_operator(cₜ₊₁, vₜ₊₁, Rₜ, Rₜ₊₁, W, πrft, pₜ, pₜ₊₁, model_params)

    # ϵ = make_ϵ(Emucₜ₊₁ , V̄ₜ₊₁ , Rₜ , W, πrft, pₜ , model_params)
    
    # shares = make_shares(Tvₜ , model_params.η, model_params.ψ)

    shares = make_shares(Tvₜ , model_params.η, model_params.ψ)

    ϵ = make_ϵ(Emucₜ₊₁ , V̄ₜ₊₁ , shares, Rₜ , W, πrft, pₜ , model_params)[1]

    return household(Kgaₜ , Kgcₜ , shares, Tvₜ , ϵ)

end

#########################################################################################
function make_dataset_NIPA(output, output_int, hh_muc, hh_muc_int, debt, G, p, T, model_params)
    # used to construct a dataset to visualize output
    backfilllength = 10

    # bigC = []
    bigT = []
    bigY = []
    bigG = []
    bigΠ = []
    bigN = []
    bigB = []
    big_poor = []
    big_rich = []
    big_middle = []

    backfill= Array{Float64}(undef, backfilllength )

    # fill!(backfill, output_int.X)

    # foo = [output[xxx].X for xxx in 1:T]

    # prepend!(foo, backfill)
    # append!(bigC, foo)

    ###############################################

    fill!(backfill, dot(output_int.X, p) )

    foo = [ dot(output[xxx].X, p) for xxx in 1:T]

    prepend!(foo, backfill)
    append!(bigY, foo)

    time = -9:1:T
    append!(bigT,time)

    ###############################################

    fill!(backfill, G)

    foo = [G for xxx in 1:T]

    prepend!(foo, backfill)
    append!(bigG, foo)

    ###############################################

    fill!(backfill, output_int.πrft)

    foo = [output[xxx].πrft for xxx in 1:T]
    
    prepend!(foo, backfill)
    append!(bigΠ, foo)

    ###############################################

    fill!(backfill, sum( output_int.Nj  ) )

    foo = [sum( output[xxx].Nj ) for xxx in 1:T]
        
    prepend!(foo, backfill)
    append!(bigN, foo)

    ###############################################

    fill!(backfill, model_params.B )

    foo = [debt[xxx] for xxx in 1:T]
            
    prepend!(foo, backfill)
    append!(bigB, foo)

    ###############################################

    fill!(backfill, hh_muc_int.poor )

    foo = [hh_muc[xxx].poor for xxx in 1:T]
                
    prepend!(foo, backfill)
    append!(big_poor, foo)

    ###############################################

    fill!(backfill, hh_muc_int.rich )

    foo = [hh_muc[xxx].rich for xxx in 1:T]
                
    prepend!(foo, backfill)
    append!(big_rich, foo)

    ###############################################

    fill!(backfill, hh_muc_int.middle )

    foo = [hh_muc[xxx].middle for xxx in 1:T]
                
    prepend!(foo, backfill)
    append!(big_middle, foo)

    ###############################################

    df = DataFrame(Y= bigY,
    Π = bigΠ,
    N = bigN,
    debt = bigB,
    G = bigG,
    poor_savings = big_poor,
    rich_savings = big_rich,
    middle_savings = big_middle,
    time = bigT,
    );

return df

end



function make_dataset(output, output_int, ppath, p_int, T, model_params)
    # used to construct a dataset to visualize output
    backfilllength = 10

    bigP = []
    bigμ = []
    bigT = []
    bigN = []
    bigS = []
    bigC = []
    good_index = []

    backfill= Array{Float64}(undef, backfilllength )

    for idxj = 1:model_params.Ngoods
        
        ###############################################
        # prices 
        fill!(backfill, p_int[idxj])

        foo = ppath[idxj,:]

        prepend!(foo, backfill)
    
        time = -9:1:T
    
        append!(bigP, foo)
        append!(bigT,time)
        append!(good_index, Int.(idxj.*ones(length(time))))

        ###############################################
        # markups

        fill!(backfill, output_int.markup[idxj])

        foo = [output[xxx].markup[idxj] for xxx in 1:T]

        prepend!(foo, backfill)
        append!(bigμ, foo)

        ###############################################
        # size

        fill!(backfill, output_int.Nj[idxj])

        foo = [output[xxx].Nj[idxj] for xxx in 1:T]

        prepend!(foo, backfill)
        append!(bigN, foo)

        ###############################################
        # size

        fill!(backfill, ( output_int.X[idxj] / sum(output_int.X ) ))

        foo = [( output[xxx].X[idxj] / sum(output[xxx].X ) ) for xxx in 1:T]

        prepend!(foo, backfill)
        append!(bigS, foo)

        ###############################################
        # size

        fill!(backfill, output_int.X[idxj] )

        foo = [( output[xxx].X[idxj] ) for xxx in 1:T]

        prepend!(foo, backfill)
        append!(bigC, foo)
    
    end

    df = DataFrame(good_index = good_index, 
               price = bigP,
               time = bigT,
               markup = bigμ,
               size = bigN,
               share = bigS,
               cj = bigC,
               );

    return df

end