struct NIPA{T}
    C::Array{T}
    AD::Array{T}
    income::Array{T}
    N::Array{T}
    Aprime::Array{T}
    LFP::Array{T}
end

##############################################################################

function get_aprime(asset_policy, work_policy, state_index, model_params)
    # grabs the choice of assets given states
    # shock is first index,
    # assets is second index,

    @unpack Woptions, agrid, statesize = model_params
    
    agrid = convert(Array{Float64}, agrid)
    
    asset_prime = Array{eltype(asset_policy)}(undef, statesize)
    fill!(asset_prime, 0.0) #need to fill given += operator below

    @views @inbounds for (foo, xxx) in enumerate(state_index)

        for wrk = 1:Woptions

            asset_prime[foo] += sum(agrid .* asset_policy[xxx[1], : ,xxx[2], wrk] * work_policy[xxx[1], xxx[2], wrk] )
                    # aprime choice given work and car choice
                    # which is all assets * prob its chosen sumed over 
        end

    end

    return asset_prime

end

##############################################################################

function get_consumption(Pces, W, τ_rev, R, asset_policy, work_policy, state_index, model_params)
 
    @unpack Na, Nshocks, Woptions, agrid, statesize, mc = model_params

    c = Array{eltype(R)}(undef, statesize)
    fill!(c, 0.0)

    exp_statevals = exp.(mc.state_values)

    @views @inbounds for (foo, xxx) in enumerate(state_index)

        for wrk = 1:Woptions

            wz = labor_income(exp_statevals[xxx[2]], W, wrk)
            # will return labor income depending upon how much working.
         
            c[foo] += sum( asset_policy[xxx[1], : ,xxx[2], wrk] * work_policy[xxx[1], xxx[2], wrk]
                        .* consumption(Pces, τ_rev, R*agrid[xxx[1]], agrid, wz)) 
        # given different aprim, all the consumption 

        end

    end

    return c

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

function get_laborincome(W, work_policy, state_index, model_params)
    # grabs the current shock in level terms
    # asset is first index,
    # shock is second index,
    # car is thrid index 
    @unpack statesize, Woptions, mc = model_params
    
    wz = Array{eltype(W)}(undef, statesize)
    fill!(wz,0.0)

    ef_units = similar(wz)
    fill!(ef_units,0.0)

    for (foo, xxx) in enumerate(state_index)

        ef_units[foo] = exp.(mc.state_values[xxx[2]])

        for wrk = 1:Woptions

            wz[foo] += labor_income(ef_units[foo], W, wrk).* work_policy[xxx[1], xxx[2], wrk]
            # will return labor income depending upon how much working.

        end

    end

    return wz, ef_units

end

##############################################################################

function get_labor(W, work_policy, state_index, distribution, model_params)

    @unpack statesize = model_params

    wz, ef_units = get_laborincome(W, work_policy, state_index, model_params)

    labor_supply = Array{eltype(W)}(undef, statesize)

    fill!(labor_supply, 0.0) #need to fill given += operator below

    for (foo, xxx) in enumerate(state_index)

        labor_supply[foo] += work_policy[xxx[1], xxx[2], 1]

    end

    LFP = sum(labor_supply .* distribution, dims = 1)[1]

    N = sum( ef_units .* labor_supply .* distribution, dims = 1)[1]

    return labor_supply, N, LFP, wz, ef_units
    # matrix that is of size Na*Nshocks*Ncars

end


##############################################################################

function get_astate(state_index, model_params)
    # grabs the current asset level
    # shock is first index,
    # assets is second index,
    # car is thrid index 
    @unpack statesize, agrid = model_params
    
    asset_state = Array{eltype(agrid)}(undef, statesize)

    for (foo, xxx) in enumerate(state_index)
        
        asset_state[foo] = agrid[xxx[1]]

    end

    return asset_state

end


##############################################################################

function aggregate(Pces, W, τ_rev, R, hh, distribution, TFP, model_params; display = false)

    # organization...

    @unpack state_index, λ = distribution

    #####
    labor_supply, N, LFP, wz, ef_units = get_labor(W, hh.work_policy, state_index, λ, model_params)

    #####
    # Get stuff from hh side

    a = get_astate(state_index, model_params)
    # assets today

    aprime = get_aprime(hh.asset_policy, hh.work_policy, state_index, model_params)
    # # assets tomorrow Noption * statesize as assets are contingent on car choice

    c = get_consumption(Pces, W, τ_rev, R, hh.asset_policy, hh.work_policy, state_index, model_params)

    Aprime = sum(aprime .* λ, dims = 1)[1]

    # aggregate asset demand 

    A = sum(a .* λ, dims = 1)[1]
    # aggregate asset positoin entering the period

    NetA = -(R * A) + Aprime

    C = sum( c .* λ, dims = 1)[1]
    #aggregate nondurable consumption
    # then multiply the prob of being in that state and sum to aggregate 

    AD = C + NetA
    # this is aggregate demand. So how many coconuts do you need which 
    # includes consumption and then stuff to save /lend out

    ####
    # then compute GDP like measure

    expenditure = Pces *( C + NetA - τ_rev)
    # to get expenditure to line up with production and income,
    # need to put everything in same units, mulitipying by P index does this
    # net off tariff revenue too

    price_production = W / TFP # comes of profit max condition that p = w / A

    production = price_production * sum( TFP .* ef_units .* labor_supply .* λ, dims = 1)[1] 

    income = sum( wz.* λ, dims = 1)[1] 
    # aggregate labor income 
    # take shock and multiply by prob of being in that state, sum to aggregate

    output_stats = NIPA(
        [C],
        [AD],
        [income],
        [N],
        [Aprime],
        [LFP]
        )

    if display == true
        digits = 6

        println("")
        println("NIPA")
        println("---------------------------------")
        println("")
        println("Aggregate Non-Durable Consumption")
        println(round(C, digits = digits))
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
        println("")
        println("Aggregate Participation Rate")
        println(round(LFP, digits = digits))

    end



    return output_stats

end

##############################################################################
##############################################################################

function make_dataset(output, output_int, Tperiods)

    backfilllength = 10
    bigT = []
    bigC = []
    bigN = []
    bigLFP = []
    country_index = []

    backfill= Array{Float64}(undef, backfilllength )

    for cnt = 1:trade_params().Ncntry

        ###############################################
        # Counsumption
        fill!(backfill, output_int[cnt].C[1])

        C = [output[cnt][xxx].C[1] for xxx in 1:Tperiods]

        prepend!(C, backfill)
        append!(bigC, C)

        ###############################################
        # Time and country index

        time = -9:1:(length(C)-10)
    

        append!(bigT,time)
        append!(country_index, Int.(cnt.*ones(length(time))))
    
        ###############################################
        # Labor Supply
        fill!(backfill, output_int[cnt].N[1])

        N = [output[cnt][xxx].N[1] for xxx in 1:Tperiods]

        prepend!(N, backfill)
        append!(bigN, N)


        ###############################################
        # LFP
        fill!(backfill, output_int[cnt].LFP[1])

        LFP = [output[cnt][xxx].LFP[1] for xxx in 1:Tperiods]

        prepend!(LFP, backfill)
        append!(bigLFP, LFP)
    end

    df = DataFrames.DataFrame(country_index = country_index, 
               consumption = bigC,
               labor_supply = bigN,
                time = bigT,
                LFP = bigLFP
               );

    return df

end