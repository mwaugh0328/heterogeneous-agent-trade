##########################################################################
##########################################################################

function trade_equilibrium(w, AD, N, τ_revenue, trade_params; output = "solver")
    # constructs zero function, takes wages, demand, tariff revenue
    # returns diffrence between procution and demand and guessed tariff transfer
    # and relized tariff transfer
    #
    # If output = "all" then returns all trade statistics

    @unpack A, Ncntry = trade_params

    p, Pindex = goods_prices(w, trade_params)

    trade = trade_flows(p, Pindex, AD, trade_params)

    value_production = similar(trade.world_demand)

    for xxx = 1:Ncntry

        value_production[xxx] = A[xxx] * N[xxx] * p[xxx,xxx]

    end

    if output == "solver"

        net_demand = value_production .- trade.world_demand
        # left hand side is production of each commodity
        # right hand side is demand of each commodity by all countries

        τ_zero = τ_revenue .- sum(trade.τ_revenue, dims = 2)[:]

        return [net_demand ; τ_zero]

    elseif output == "all"

        return trade

    end

end

##########################################################################
##########################################################################



