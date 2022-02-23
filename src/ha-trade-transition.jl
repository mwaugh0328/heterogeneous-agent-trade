function transition_path_FixedR(x, FixedR, Rint, τ_path, dist, hh, model_params; display = false)
    # for fixed R...still uses basic functions...would like a nicer multiple dispatch setup
    # here the asset market does not clear
    
    Ncntry = τ_path[1].Ncntry
        
    Tperiods = size(τ_path)[1]

    Wvec = x[ 1: Ncntry*Tperiods]
    τ_rev_vec = x[Ncntry*Tperiods + 1: 2*Ncntry*Tperiods]

    Rvec = Array{eltype(FixedR)}(undef, (Ncntry*Tperiods - 1))
    fill!(Rvec, FixedR)

    Rvec = vcat(Rint, Rvec)

    if display == false

        goods_demand, asset_demand, output_stats = transition_path(Wvec, τ_rev_vec, Rvec, τ_path, dist, hh, model_params)

        return goods_demand

    else

        goods_demand, asset_demand, output_stats = transition_path(Wvec, τ_rev_vec, Rvec, τ_path, dist, hh, model_params)

        return goods_demand, asset_demand, output_stats

    end

end


function transition_path(x, Rint, τ_path, dist, hh, model_params; display = false)
    #using mulitple dispatch here
    
    Ncntry = τ_path[1].Ncntry
        
    Tperiods = size(τ_path)[1]

    Wvec = x[ 1: Ncntry*Tperiods]
    τ_rev_vec = x[Ncntry*Tperiods + 1: 2*Ncntry*Tperiods]
    
    
    Rvec = x[2*Ncntry*Tperiods + 1 : end]
    Rvec = vcat(Rint, Rvec)

    if display == false

        goods_demand, asset_demand, output_stats = transition_path(Wvec, τ_rev_vec, Rvec, τ_path, dist, hh, model_params)

        return vcat(goods_demand, asset_demand)

    else

        goods_demand, asset_demand, output_stats = transition_path(Wvec, τ_rev_vec, Rvec, τ_path, dist, hh, model_params)

        return vcat(goods_demand, asset_demand), output_stats

    end

end

#####################################################################################################
#####################################################################################################

function transition_path(Wvec, τ_rev_vec, Rvec, τ_path, dist, hh, model_params)
    # main version to be used in solving for eq.

    @assert length(Wvec) == length(τ_rev_vec)

    Ncntry = τ_path[1].Ncntry
        
    Tperiods = size(τ_path)[1]
    # \tau_path is a structure of the trade params. should have every fundemental
    # there for each period.

    # then this just organizes things so stuff will
    # be of size Ncntry, Tperiods

    W = reshape(Wvec, Ncntry, Tperiods)
    τ_rev = reshape(τ_rev_vec, Ncntry, Tperiods)

    Pces = similar(W)
    AD = similar(W)
    N = similar(W)
    AP = similar(W)

    # this stuff below is for finacial integration
    if length(Rvec) == Ncntry*Tperiods

        R = reshape(Rvec, Ncntry, Tperiods)

    else

        @assert length(Rvec) == Tperiods

        R = reshape(repeat(Rvec, inner = Ncntry, outer = 1), Ncntry, Tperiods)

        @assert size(R) == size(W)

    end

    output_stats = Array{Array{NIPA{eltype(W)}}}(undef, Ncntry)

    net_demand = Array{eltype(Wvec)}(undef, 2*Ncntry, Tperiods)
    # The output vectov, net demand is 2*Ncountry bc we are 
    # checking if output and tariff revenue is consistent.

    for fwdate = 1:Tperiods
        # create prices for every period

        Pces[:, fwdate] = goods_prices(W[:, fwdate], τ_path[fwdate])[2]

    end

    Threads.@threads for nct = 1:Ncntry 
        # pass through the backward/forward step.
        # this is an effective place to multithread (if not distribute)

        TFP = [τ_path[xxx].A[nct] for xxx = 1:Tperiods]

        AD[nct, :], N[nct, :], AP[nct,:], output_stats[nct] = foward_backward(Pces[nct,:], W[nct,:], 
                                    τ_rev[nct,:], R[nct,:], dist[nct], hh[nct].Tv, TFP, model_params) 

    end

    AP = AP[:,1:end-1]

    # Then given how each country behaves, we can bring it together through the 
    # trade eq....

    for fwdate = 1:Tperiods

        net_demand[:, fwdate] = trade_equilibrium(W[:, fwdate], AD[:, fwdate], 
                N[:, fwdate], τ_rev[:, fwdate], τ_path[fwdate])

    end

    net_demand = net_demand[:]

    # Again this is the situation with the
    # finacial autarky or global integration
    if length(Rvec) == Ncntry*Tperiods

        AP = AP[:]
    
    else
        # sum across the country dimension
        AP = sum(AP, dims = 1)

        AP = AP[:]

    end

    # this is the finacial autarky situation. 
    return net_demand, AP , output_stats

end

#####################################################################################################
#####################################################################################################

function foward_backward(Pces, W, τ_rev, R, dist, Tvend, TFP, model_params)
    # Foward backward approach to compute path.

    @unpack Na, Nshocks, Woptions, mc, β, σa, σw = model_params

    Tperiods = length(W)

    hh = Array{household{eltype(W)}}(undef, Tperiods)
    # A structure of households (policies and Vs) for each time period

    output_stats = Array{NIPA{eltype(W)}}(undef, Tperiods)
    # structure of statistics, for each time period.

    AD = Array{eltype(W)}(undef, Tperiods)
    N = similar(AD)

    u = Array{eltype(W)}(undef, Na, Na, Nshocks, Woptions)

    #####################################################################################################
    # This is the backward step: sovle hh problem at Tperiods then use bellman operator to work backwards

    for bwdate = Tperiods:-1:1
        # R is fixed, but it's value depends upon price level last period, that's the value of stuff
        # invested. So you need to take it into account.

        make_utility!(u, Pces[bwdate], W[bwdate], τ_rev[bwdate], R[bwdate], model_params)

        if bwdate == Tperiods

            hh[bwdate] = bellman_operator_policy(Tvend, u , mc.p, β, σa, σw);

        else

            hh[bwdate] = bellman_operator_policy(hh[bwdate + 1].Tv, u , mc.p, β, σa, σw);
            # + 1 in hh because we are working backwards

        end
        
        # use the relationship where
        # Tv_{date} = max [u_{date} + BETV_{date + 1} ]
        # and policy functions for date T are the arg max.

    end

    #####################################################################################################
    # This is the foward step: Start from the inital distribution λ, make the Q given policies hh at that 
    # date, then compute λ_{t+1}

    @unpack Q, λ, state_index = dist
    
    @views for fwdate = 1:Tperiods

        output_stats[fwdate] = aggregate(Pces[fwdate], W[fwdate], τ_rev[fwdate], R[fwdate],
                                     hh[fwdate], distribution(Q, λ, state_index), TFP[fwdate], model_params)

        λ = push_foward(λ, Q, state_index, hh[fwdate] , model_params)

    end

    AD = [output_stats[xxx].AD[1] for xxx in 1:Tperiods]
    N  = [output_stats[xxx].N[1] for xxx in 1:Tperiods]
    AP = [output_stats[xxx].Aprime[1] for xxx in 1:Tperiods]

    return AD, N, AP, output_stats

end

#####################################################################################################
#####################################################################################################

function push_foward(λ, Q, state_index, household, model_params)
# Pushes the economy forward
# (1) Take policy functions and a distribution L -> aggregates today.
# Policy functions -> Transition probability Q
# Distribution today + Q -> Distribution tomorrow. 

# make transition matrix

make_Q!(Q, state_index, household.asset_policy, household.work_policy, model_params)

# Then push the distribution forward
return law_of_motion(λ , transpose(Q))

end

#####################################################################################################
#####################################################################################################

function collect_end_conditions(W::Array{T}, τ_rev::Array{T}, 
                            R::Array{T}, model_params, trade_params) where T
    # grabs end conditions

    @unpack Ncntry = trade_params

    Pces = goods_prices(W, trade_params)[2]
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    hh = Array{household{eltype(W)}}(undef, Ncntry)

    Threads.@threads for ncnty = 1:Ncntry  # do this for each country.
                            # this is the place to use distributed. hh problems can be solved independtly
        hh[ncnty] = solve_household_problem(Pces[ncnty], W[ncnty], τ_rev[ncnty], R[ncnty],  model_params)
        
    end

    return hh

end

function collect_end_conditions(W::Array{T}, τ_rev::Array{T}, R::T, model_params, trade_params) where T
    # grabs end conditions

    @unpack Ncntry = trade_params

    Pces = goods_prices(W, trade_params)[2]
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    hh = Array{household{eltype(W)}}(undef, Ncntry)

    Threads.@threads for ncnty = 1:Ncntry  # do this for each country.
                            # this is the place to use distributed. hh problems can be solved independtly
        hh[ncnty] = solve_household_problem(Pces[ncnty], W[ncnty], τ_rev[ncnty], R,  model_params)
        
    end

    return hh

end

##########################################################################
##########################################################################


function collect_intial_conditions(W::Array{T}, τ_revenue::Array{T}, 
                                R::Array{T}, model_params, trade_params) where T
    # grabs initial conditions for 
    # to compute 

    @unpack Na, statesize = model_params
    @unpack Ncntry = trade_params

    Pces = goods_prices(W, trade_params)[2]
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    dist = Array{distribution{eltype(W)}}(undef, Ncntry)

    Threads.@threads for ncnty = 1:Ncntry  # do this for each country.

        dist[ncnty] = compute_eq(Pces[ncnty], W[ncnty], τ_revenue[ncnty], R[ncnty],  model_params)[2];

    end

    return dist

end

function collect_intial_conditions(W::Array{T}, τ_revenue::Array{T},
                         R::T, model_params, trade_params) where T
    # grabs initial conditions for 
    # to compute 

    @unpack Na, statesize = model_params
    @unpack Ncntry = trade_params

    Pces = goods_prices(W, trade_params)[2]
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    dist = Array{distribution{eltype(W)}}(undef, Ncntry)

    Threads.@threads for ncnty = 1:Ncntry  # do this for each country.

        dist[ncnty] = compute_eq(Pces[ncnty], W[ncnty], τ_revenue[ncnty], R,  model_params)[2];

    end

    return dist

end


