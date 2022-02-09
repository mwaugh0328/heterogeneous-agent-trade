function ha_trade_equilibrium(W, τ_revenue, R, model_params, trade_params; output = "solver")
    # Combines the demand side with the trade side.

    @unpack Ncntry = trade_params

    AD = Array{eltype(W)}(undef, Ncntry)
    Pces = similar(AD)
    N = similar(AD)
    LFP = similar(AD)


    foo, Pces = goods_prices(W, trade_params)
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    for xxx = 1:Ncntry  # do this for each country.
                        # this is the place to use distributed. hh problems can be solved independtly

        hh, dist = compute_eq(Pces[xxx], W[xxx], τ_revenue[xxx], R[xxx],  model_params);
        # This is the BHA part. Given prices, solve hh problem and construct stationary distribution L
        # internaltionally all the matters from the hh perspective is the price index.

        output_stats[xxx] = aggregate(Pces[xxx], W[xxx], τ_revenue[xxx], R, trade_params.A[xxx], 
                model_params, ap, cp, wp, state_index, L; display = true);

    end
    
    # We have each countries total demand of goods: AD
    # And each countries total supply of goods: A*N
    
    # Pass supply and demand to the trade block and see if they match up.

    if output == "solver"

        return compute_trade_eq(W, AD, N, τ_revenue, trade_params)

    elseif output == "all"

        return compute_trade_eq(W, AD, N, τ_revenue, trade_params), output_stats
    end

end