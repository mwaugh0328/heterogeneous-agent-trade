##########################################################################
##########################################################################

function trade_equilibrium(w, AD, N, τrev, trade_params; output = "solver")
    # constructs zero function, takes wages, demand, tariff revenue
    # returns diffrence between procution and demand and guessed tariff transfer
    # and relized tariff transfer
    #
    # If output = "all" then returns all trade statistics

    @unpack A, Ncntry = trade_params

    p, Pindex = goods_prices(w, trade_params)

    trade = trade_flows(p, Pindex, AD, trade_params)

    value_production = Array{eltype(w)}(undef, Ncntry)
    net_demand = similar(value_production)
    τ_zero = similar(value_production)

    for xxx = 1:Ncntry

        value_production[xxx] = A[xxx] * N[xxx] * p[xxx,xxx]

    end

    @unpack world_demand, τ_revenue = trade

     net_demand .= value_production .- world_demand
        # left hand side is production of each commodity
        # right hand side is demand of each commodity by all countries

    τ_zero .= τrev .- sum(τ_revenue, dims = 2)[:]

    return vcat(net_demand , τ_zero)

end

##########################################################################
##########################################################################



