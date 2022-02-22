# Note how multiple dispatch is being used here...

##########################################################################
##########################################################################

function ha_trade_equilibrium(x, model_params, trade_params; display = false)
    #again using mulitple dispatch here...size of x determines
    
    @unpack Ncntry = trade_params

    if length(x) == 3*Ncntry #Financial Autarky case

        goods_demand , asset_demand, output_stats = ha_trade_equilibrium(x[1:Ncntry], x[Ncntry+1 : 2*Ncntry], 
                    x[2*Ncntry + 1 : end], model_params, trade_params)

    elseif length(x) == 2*Ncntry + 1 #Financial Integration case

        goods_demand , asset_demand, output_stats = ha_trade_equilibrium(x[1:Ncntry], x[Ncntry+1 : 2*Ncntry], 
                    x[2*Ncntry + 1], model_params, trade_params)

    end

    if display == false

        return vcat(goods_demand , asset_demand)

    else

        return vcat(goods_demand , asset_demand ), output_stats

    end

end

##########################################################################
##########################################################################

function ha_trade_equilibrium(x, R, model_params, trade_params; display = false)
    #again using mulitple dispatch here for fixed R
    
    @unpack Ncntry = trade_params

    goods_demand, asset_demand, output_stats = ha_trade_equilibrium(x[1:Ncntry], x[Ncntry+1 : 2*Ncntry], 
                    R, model_params, trade_params)


    if display == false

        return goods_demand

    else

        return goods_demand, output_stats

    end

end


##########################################################################
##########################################################################

function ha_trade_equilibrium(W::Array{T}, τ_revenue::Array{T}, R::Array{T}, model_params, trade_params) where T
    # Combines the demand side with the trade side.
    #
    # Using multiple dispatch here to determine the level of finacial integration.
    # If type of R is array -> each country has it's own rate -> in fincial autarky, so the bond market
    # needs to clear within the country

    @unpack Ncntry, A = trade_params

    output_stats = Array{NIPA{eltype(W)}}(undef, Ncntry)
    hh = Array{household{eltype(W)}}(undef, Ncntry)
    dist = Array{distribution{eltype(W)}}(undef, Ncntry)

    Pces = goods_prices(W, trade_params)[2]

    
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    Threads.@threads for xxx = 1:Ncntry  # do this for each country.
                        # this is the place to use distributed or pmap
                        # hh problems can be solved independtly

        hh[xxx], dist[xxx] = compute_eq(Pces[xxx], W[xxx], τ_revenue[xxx], R[xxx],  model_params);
        # This is the BHA part. Given prices, solve hh problem and construct stationary distribution L
        # internaltionally all the matters from the hh perspective is the price index.

        output_stats[xxx] = aggregate(Pces[xxx], W[xxx], τ_revenue[xxx], R[xxx],
                                     hh[xxx], dist[xxx], A[xxx], model_params)

    end
    
    # We have each countries total demand of goods: AD
    # And each countries total supply of goods: A*N
    
    # Pass supply and demand to the trade block and see if they match up.

    AD = [output_stats[cnt].AD[1] for cnt in 1:Ncntry]

    Nf = [output_stats[cnt].N[1] for cnt in 1:Ncntry]

    asset_net_demand = [output_stats[cnt].Aprime[1] for cnt in 1:Ncntry]
    # then this is the asset market clearing condition

    trade_net_demand = trade_equilibrium(W, AD, Nf, τ_revenue, trade_params)

    return trade_net_demand , asset_net_demand, output_stats

end

##########################################################################
##########################################################################

function ha_trade_equilibrium(W::Array{T}, τ_revenue::Array{T}, R::T, model_params, trade_params) where T
    # Combines the demand side with the trade side.
    #
    # Using multiple dispatch here to determine the level of finacial integration.
    # If type of R is singlton of type T -> there is world interest rate for which the bond market
    # needs to clear at that one rate. 

    # I think not hard typing, using T helps if autodiff is employed

    @unpack Ncntry, A = trade_params

    output_stats = Array{NIPA{eltype(W)}}(undef, Ncntry)
    hh = Array{household{eltype(W)}}(undef, Ncntry)
    dist = Array{distribution{eltype(W)}}(undef, Ncntry)

    Pces = goods_prices(W, trade_params)[2]

    
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    Threads.@threads for xxx = 1:Ncntry  # do this for each country.
                        # this is the place to use distributed or pmap
                        # hh problems can be solved independtly

        hh[xxx], dist[xxx] = compute_eq(Pces[xxx], W[xxx], τ_revenue[xxx], R,  model_params);
        # This is the BHA part. Given prices, solve hh problem and construct stationary distribution L
        # internaltionally all the matters from the hh perspective is the price index.

        output_stats[xxx] = aggregate(Pces[xxx], W[xxx], τ_revenue[xxx], R,
                        hh[xxx], dist[xxx], A[xxx], model_params)

    end
    
    # We have each countries total demand of goods: AD
    # And each countries total supply of goods: A*N
    
    # Pass supply and demand to the trade block and see if they match up.

    AD = [output_stats[cnt].AD[1] for cnt in 1:Ncntry]

    Nf = [output_stats[cnt].N[1] for cnt in 1:Ncntry]

    asset_net_demand = sum([output_stats[cnt].Aprime[1] for cnt in 1:Ncntry])
    # then this is the asset market clearing condition
    # sum across world asset demand---given how mulitiple dispatch is working

    trade_net_demand = trade_equilibrium(W, AD, Nf, τ_revenue, trade_params)

    return trade_net_demand , asset_net_demand , output_stats

end

##########################################################################
##########################################################################

function ha_trade_equilibrium_pmap(x, model_params, trade_params)
    #again using mulitple dispatch here...size of x determines
    
    @unpack Ncntry = trade_params

    if length(x) == 3*Ncntry #Financial Autarky case

        outvec = ha_trade_equilibrium_pmap(x[1:Ncntry], x[Ncntry+1 : 2*Ncntry], x[2*Ncntry + 1 : end], model_params, trade_params)

    elseif length(x) == 2*Ncntry + 1 #Financial Integration case

        outvec = ha_trade_equilibrium_pmap(x[1:Ncntry], x[Ncntry+1 : 2*Ncntry], x[2*Ncntry + 1], model_params, trade_params)

    end

    return outvec

end

##########################################################################

function ha_trade_equilibrium_pmap(W::Array{T}, τ_revenue::Array{T}, R::Array{T}, model_params, trade_params) where T
    # Combines the demand side with the trade side.
    #
    # Using multiple dispatch here to determine the level of finacial integration.
    # If type of R is array -> each country has it's own rate -> in fincial autarky, so the bond market
    # needs to clear within the country

    @unpack Ncntry, A = trade_params

    output_stats = Array{NIPA{eltype(W)}}(undef, Ncntry)

    Pces = goods_prices(W, trade_params)[2]

    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    x = Array{Tuple{eltype(W), eltype(W), eltype(W), eltype(W) }}(undef, Ncntry)

    for xxx = 1:Ncntry

        x[xxx] = (Pces[xxx], W[xxx], τ_revenue[xxx], R[xxx])

    end

    g(x) = compute_eq(x[1], x[2], x[3], x[4], model_params)

    outvec = pmap(g,x)

    for xxx = 1:Ncntry  # do this for each country.
                        # this is the place to use distributed or pmap
                        # hh problems can be solved independtly

        # hh[xxx], dist[xxx] = compute_eq(Pces[xxx], W[xxx], τ_revenue[xxx], R[xxx],  model_params);
        # # This is the BHA part. Given prices, solve hh problem and construct stationary distribution L
        # # internaltionally all the matters from the hh perspective is the price index.

        output_stats[xxx] = aggregate(Pces[xxx], W[xxx], τ_revenue[xxx], R[xxx],
                    outvec[xxx][1], outvec[xxx][2], A[xxx], model_params)

    end
    
    # We have each countries total demand of goods: AD
    # And each countries total supply of goods: A*N
    
    # Pass supply and demand to the trade block and see if they match up.

    AD = [output_stats[cnt].AD[1] for cnt in 1:Ncntry]

    Nf = [output_stats[cnt].N[1] for cnt in 1:Ncntry]

    asset_net_demand = [output_stats[cnt].Aprime[1] for cnt in 1:Ncntry]
    # then this is the asset market clearing condition

    trade_net_demand = trade_equilibrium(W, AD, Nf, τ_revenue, trade_params)

    return vcat(trade_net_demand , asset_net_demand )

end

##########################################################################
##########################################################################


function ha_trade_equilibrium_pmap(W::Array{T}, τ_revenue::Array{T}, R::T, model_params, trade_params) where T
    # Combines the demand side with the trade side.
    #
    # Using multiple dispatch here to determine the level of finacial integration.
    # If type of R is array -> each country has it's own rate -> in fincial autarky, so the bond market
    # needs to clear within the country

    @unpack Ncntry, A = trade_params

    output_stats = Array{NIPA{eltype(W)}}(undef, Ncntry)

    Pces = goods_prices(W, trade_params)[2]

    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    x = Array{Tuple{eltype(W), eltype(W), eltype(W)}}(undef, Ncntry)

    for xxx = 1:Ncntry

        x[xxx] = (Pces[xxx], W[xxx], τ_revenue[xxx])

    end

    g(x) = compute_eq(x[1], x[2], x[3], R, model_params)

    outvec = pmap(g,x)

    for xxx = 1:Ncntry  # do this for each country.
                        # this is the place to use distributed or pmap
                        # hh problems can be solved independtly

        # hh[xxx], dist[xxx] = compute_eq(Pces[xxx], W[xxx], τ_revenue[xxx], R[xxx],  model_params);
        # # This is the BHA part. Given prices, solve hh problem and construct stationary distribution L
        # # internaltionally all the matters from the hh perspective is the price index.

        output_stats[xxx] = aggregate(Pces[xxx], W[xxx], τ_revenue[xxx], R,
                    outvec[xxx][1], outvec[xxx][2], A[xxx], model_params)

    end
    
    # We have each countries total demand of goods: AD
    # And each countries total supply of goods: A*N
    
    # Pass supply and demand to the trade block and see if they match up.

    AD = [output_stats[cnt].AD[1] for cnt in 1:Ncntry]

    Nf = [output_stats[cnt].N[1] for cnt in 1:Ncntry]

    asset_net_demand = sum([output_stats[cnt].Aprime[1] for cnt in 1:Ncntry])
    # then this is the asset market clearing condition
    # sum across world asset demand---given how mulitiple dispatch is working

    trade_net_demand = trade_equilibrium(W, AD, Nf, τ_revenue, trade_params)

    return vcat(trade_net_demand , asset_net_demand )

end

##########################################################################

function unpack_solution(solution, Ncntry)

    W = solution[1:Ncntry]

    τ_rev =  solution[Ncntry+1 : 2*Ncntry]
    
    R = solution[2*Ncntry + 1 : end]

    return W, τ_rev, R
    
end

function unpack_solution(solution, R, Ncntry)

    W = solution[1:Ncntry]

    τ_rev =  solution[Ncntry+1 : 2*Ncntry]
    
    return W, τ_rev, R
    
end