struct NIPA
    PC::Float64 #consumption
    M::Float64 #imports
    X::Float64 #exports
    income::Float64
    production::Float64
    N::Float64
    Aprime::Float64
    NetA::Float64
end

struct trade
    bilateral_imports::Array{Float64} # asset_policy
    bilateral_πprob::Array{Float64} # choice probabilities
end

##############################################################################

function get_aprime(asset_policy, πprob, state_index, model_params)
    # grabs the choice of assets given states

    @unpack Na, Nshocks, Ncntry = model_params
    
    asset_prime = Array{eltype(asset_policy)}(undef, Na*Nshocks)
    fill!(asset_prime, 0.0) #need to fill given += operator below

    for (foo, xxx) in enumerate(state_index)

        for cntry = 1:Ncntry

            asset_prime[foo] += asset_policy[xxx[1], xxx[2], cntry]*πprob[xxx[1], xxx[2], cntry]

        end

    end

    return asset_prime

end

##############################################################################

function get_consumption(R, W, asset_policy, πprob, state_index, model_params)
    # grabs the choice of assets given states

    @unpack Na, Nshocks, Ncntry, mc, agrid = model_params
    
    pconsumption = Array{eltype(asset_policy)}(undef, Na*Nshocks)
    fill!(pconsumption, 0.0) #need to fill given += operator below

    shocks = exp.(mc.state_values)

    for (foo, xxx) in enumerate(state_index)

        for cntry = 1:Ncntry

            pconsumption[foo] +=  ( -asset_policy[xxx[1], xxx[2], cntry] 
                    + R*agrid[xxx[1]] + W*shocks[xxx[2]] ) * πprob[xxx[1], xxx[2], cntry]

        end
        
    end

    return pconsumption

end

##############################################################################

function get_trade(R, W, asset_policy, πprob, state_index, model_params)
    # grabs the choice of assets given states

    @unpack Na, Nshocks, Ncntry, mc, agrid = model_params
    
    c_by_variety = Array{eltype(asset_policy)}(undef, Na*Nshocks, Ncntry)
    variety_share = Array{eltype(asset_policy)}(undef, Na*Nshocks, Ncntry)

    fill!(c_by_variety, 0.0) #need to fill given += operator below
    fill!(variety_share, 0.0) #need to fill given += operator below

    shocks = exp.(mc.state_values)

    for (foo, xxx) in enumerate(state_index)

        for cntry = 1:Ncntry

            c_by_variety[foo, cntry] +=  ( -asset_policy[xxx[1], xxx[2], cntry] 
                    + R*agrid[xxx[1]] + W*shocks[xxx[2]] ) * πprob[xxx[1], xxx[2], cntry]

            variety_share[foo, cntry] += πprob[xxx[1], xxx[2], cntry]

        end
        
    end

    return c_by_variety, variety_share

end

##############################################################################

function get_distribution(state_index, stationary_distribution)
    # returns distribution of assets across grid

    asset_state = [xxx[1] for xxx in state_index]
        # asset should be second index

    asset_distribution = [sum(stationary_distribution[asset_state .== xxx]) for xxx in unique(asset_state)]

    return asset_distribution

end

##############################################################################

function get_laborincome(W, state_index, model_params)

    @unpack Na, Nshocks, mc = model_params
    
    wz = Array{eltype(W)}(undef, Na*Nshocks)
    fill!(wz,0.0)

    ef_units = similar(wz)
    fill!(ef_units,0.0)

    for (foo, xxx) in enumerate(state_index)

        ef_units[foo] = exp.(mc.state_values[xxx[2]])

        wz[foo] += labor_income(ef_units[foo], W)

    end

    return wz, ef_units

end

##############################################################################

function get_astate(state_index, model_params)

    @unpack Na, Nshocks, agrid = model_params
    
    asset_state = Array{eltype(agrid)}(undef, Na*Nshocks)

    for (foo, xxx) in enumerate(state_index)
        
        asset_state[foo] = agrid[xxx[1]]

    end

    return asset_state

end


##############################################################################

function aggregate(R, W, p, country, household, distribution, model_params; display = false)

    # organization...

    @unpack state_index, λ = distribution
    @unpack πprob, asset_policy = household
    @unpack TFP = model_params

    #####
    wz, ef_units = get_laborincome(W, state_index, model_params)

    N = dot(ef_units, λ)

    #####
    # Get stuff from hh side

    a = get_astate(state_index, model_params)
    # assets today

    aprime = get_aprime(asset_policy, πprob, state_index, model_params)
    # # assets tomorrow Noption * statesize as assets are contingent on car choice

    pc = get_consumption(R, W, asset_policy, πprob, state_index, model_params)

    c_by_variety, variety_share = get_trade(R, W, asset_policy, πprob, state_index, model_params)

    Aprime = dot(aprime, λ)

    # aggregate asset demand 

    A = dot(a , λ)
    # aggregate asset positoin entering the period

    NetA = -(R * A) + Aprime

    PC = dot(pc , λ)
    #aggregate consumption

    imports = sum(c_by_variety[:, 1:end .!= country], dims = 2)
    # exclude home country imports

    M = dot(imports , λ)

    bilateral_imports = sum( c_by_variety .* λ, dims = 1)
    bilateral_πprob = sum( variety_share  .* λ, dims = 1)

    production = p[country] * TFP[country]* N

    home_consumption = sum(c_by_variety[:, 1:end .== country], dims = 2)

    X = production - dot(home_consumption , λ)

    ####
    # then compute GDP like measure

    income = dot(wz, λ)

    expenditure = PC + X - M
    # aggregate labor income 
    # is this correct price

    if display == true
        digits = 6

        println("")
        println("NIPA")
        println("---------------------------------")
        println("")
        println("Aggregate Consumption")
        println(round(PC, digits = digits))
        println("")
        println("Aggregate Imports")
        println(round(M, digits = digits))
        println("")
        println("Aggregate Exports")
        println(round(X, digits = digits))
        println("")
        println("Aggregate Net Asset Position")
        println(round(NetA, digits = digits))
        println("")
        println("GDP, Expenditure Side")
        println(round(expenditure, digits = digits))
        println("")
        println("GDP, Income Side")
        println(round(income, digits = digits))
        println("")
        println("GDP, Produciton Side")
        println(round(production, digits = digits))

    end

    return NIPA(PC, M, X, income, production, N, Aprime, NetA), trade(bilateral_imports, bilateral_πprob)

end

##############################################################################
##############################################################################
