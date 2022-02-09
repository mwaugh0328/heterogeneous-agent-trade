function ha_trade_equilibrium(W, τ_revenue, R, model_params, trade_params)
    # Combines the demand side with the trade side.

    @unpack Ncntry, A = trade_params

    output_stats = Array{NIPA{eltype(W)}}(undef, Ncntry)
    hh = household{eltype(W)}
    dist = distribution{eltype(W)}

    Pces = goods_prices(W, trade_params)[2]

    
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    for xxx = 1:Ncntry  # do this for each country.
                        # this is the place to use distributed or pmap
                        # hh problems can be solved independtly

        hh, dist = compute_eq(Pces[xxx], W[xxx], τ_revenue[xxx], R[xxx],  model_params);
        # This is the BHA part. Given prices, solve hh problem and construct stationary distribution L
        # internaltionally all the matters from the hh perspective is the price index.

        output_stats[xxx] = aggregate(Pces[xxx], W[xxx], τ_revenue[xxx], R[xxx],
                                     hh, dist, A[xxx], params)

    end
    
    # We have each countries total demand of goods: AD
    # And each countries total supply of goods: A*N
    
    # Pass supply and demand to the trade block and see if they match up.

    AD = [output_stats[cnt].AD[1] for cnt in 1:Ncntry]

    N = [output_stats[cnt].N[1] for cnt in 1:Ncntry]

    asset_net_demand = [output_stats[cnt].Aprime[1] for cnt in 1:Ncntry]
    # then this is the asset market clearing condition

    trade_net_demand = trade_equilibrium(W, AD, N, τ_revenue, trade_params)

    return [ trade_net_demand ; asset_net_demand ]

end