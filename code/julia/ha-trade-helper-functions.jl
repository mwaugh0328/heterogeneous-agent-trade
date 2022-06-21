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
    
    pc_by_state = Array{eltype(asset_policy)}(undef, Na*Nshocks, Ncntry)
    pcπ_by_state = Array{eltype(asset_policy)}(undef, Na*Nshocks, Ncntry)
    πprob_by_state = Array{eltype(asset_policy)}(undef, Na*Nshocks, Ncntry)

    fill!(pc_by_state , 0.0) #need to fill given += operator below
    fill!(pcπ_by_state , 0.0) #need to fill given += operator below
    fill!(πprob_by_state, 0.0) #need to fill given += operator below

    shocks = exp.(mc.state_values)

    for (foo, xxx) in enumerate(state_index)

        for cntry = 1:Ncntry

            pc = ( -asset_policy[xxx[1], xxx[2], cntry] + R*agrid[xxx[1]] + W*shocks[xxx[2]] )
            #should watch this line

            pc_by_state[foo, cntry] += pc

            pcπ_by_state[foo, cntry] +=  pc * πprob[xxx[1], xxx[2], cntry]
                    # this is p*c*π_{ij} by state

            πprob_by_state[foo, cntry] += πprob[xxx[1], xxx[2], cntry]
            # this is π_{ij} by state

        end
        
    end

    return pcπ_by_state , πprob_by_state, pc_by_state

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
    @unpack TFP, L = model_params

    #####
    wz, ef_units = get_laborincome(W, state_index, model_params)

    N = L[country] .* dot(ef_units, λ) #number of effeciency units (16)

    #####
    # Get stuff from hh side

    a = get_astate(state_index, model_params)
    # assets today

    aprime = get_aprime(asset_policy, πprob, state_index, model_params)
    # # assets tomorrow Noption * statesize as assets are contingent on car choice

    pc = get_consumption(R, W, asset_policy, πprob, state_index, model_params)

    pcπ_by_state, π_by_state = get_trade(R, W, asset_policy, πprob, state_index, model_params)[1:2]

    Aprime = L[country] .* dot(aprime, λ) 
    # asset holdings next period (17)

    # aggregate asset demand 

    A = L[country] .* dot(a , λ)
    # aggregate asset positoin entering the period

    NetA = -(R * A) + Aprime

    PC = L[country] .* dot(pc , λ)
    #aggregate consumption

    ##############
    # the Trade Stuff

    imports = sum(pcπ_by_state[:, 1:end .!= country], dims = 2)
    # exclude home country imports

    M = L[country] .* dot(imports , λ)
    # aggregate imports by value

    bilateral_imports = sum( L[country] .* pcπ_by_state.* λ, dims = 1)

    bilateral_πprob = sum( π_by_state.* λ, dims = 1)

    production = p[country] * TFP[country]* N 

    home_consumption = sum(pcπ_by_state[:, 1:end .== country], dims = 2)

    X = production - L[country] .* dot(home_consumption , λ)

    #############
    # then compute GDP like measure

    income = L[country] .* dot(wz, λ)

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
