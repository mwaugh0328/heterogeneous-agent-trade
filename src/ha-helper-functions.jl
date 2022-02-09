struct NIPA{T}
    C::Array{T}
    AD::Array{T}
    income::Array{T}
    N::Array{T}
    Aprime::Array{T}
    LFP::Array{T}
end

##############################################################################

function get_aprime(asset_policy, state_index, model_params)
    # grabs the choice of assets given states and car choice
    # shock is first index,
    # assets is second index,
    # car is thrid index 
    # but in the asset_policy it is setup with row = asset, shock is column, then thrid
    # dimension is car state, fourth dimension is if you buy or not. 

    @unpack Na, Nshocks, agrid, statesize = model_params
       
    asset_prime = Array{eltype(asset_policy)}(undef, statesize)
    fill!(asset_prime, 0.0) #need to fill given += operator below

    for (foo, xxx) in enumerate(state_index)

        asset_prime[foo] += sum(agrid .* asset_policy[xxx[1], : ,xxx[2]])
                    # aprime choice given work and car choice
                    # which is all assets * prob its chosen sumed over 
    end

    return asset_prime

end

##############################################################################

function get_consumption(Pces, W, τ_rev, R, asset_policy, model_params, state_index)
 
    @unpack Na, Nshocks, agrid, statesize, mc = model_params

    c = Array{eltype(R)}(undef, statesize)
    fill!(c, 0.0)

    for (foo, xxx) in enumerate(state_index)

        wz = labor_income(exp.(mc.state_values[xxx[2]]), W)
            # will return labor income depending upon how much working.
         
        c[foo] += sum( asset_policy[xxx[1], : ,xxx[2]] .* consumption(Pces, τ_rev, R*agrid[xxx[1]], agrid, wz)) 
        # given different aprim, all the consumption 

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

function get_laborincome(W, state_index, model_params)
    # grabs the current shock in level terms
    # asset is first index,
    # shock is second index,
    # car is thrid index 
    @unpack statesize, mc = model_params
    
    wz = Array{eltype(mc.state_values)}(undef, statesize)
    fill!(wz,0.0)

    ef_units = similar(wz)

    for (foo, xxx) in enumerate(state_index)

        ef_units[foo] = exp.(mc.state_values[xxx[2]])

        wz[foo] = labor_income(ef_units[foo], W)
            # will return labor income depending upon how much working.

    end

    return wz, ef_units

end

##############################################################################

function get_laborsupply(W, work_policy, state_index, distribution, model_params)

    @unpack Na,  Nshocks, agrid, statesize = model_params
       
    wz, ef_units = get_laborincome(W, state_index, model_params)

    labor_supply = Array{eltype(work_policy)}(undef, statesize)

    fill!(labor_supply, 0.0) #need to fill given += operator below

    for (foo, xxx) in enumerate(state_index)

        labor_supply[foo] += work_policy[xxx[1], xxx[2]]

    end

    LFP = sum(labor_supply .* distribution, dims = 1)[1]

    N = sum( ef_units .* labor_supply .* distribution, dims = 1)[1]

    return labor_supply, N, LFP
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

    @unpack state_index, L = distribution

    @unpack asset_policy = hh

    # it's setup to accomedate labor supply, but for it's not here...
    work_policy = Array{eltype(W)}(undef, model_params.Na, model_params.Nshocks)

    fill!(work_policy, 1.0)

    #####
    # Get stuff from labor side

    wz, ef_units = get_laborincome(W, state_index, model_params)

    labor_supply, N, LFP = get_laborsupply(W, work_policy, state_index, L, model_params)

    #####
    # Get stuff from hh side

    a = get_astate(state_index, model_params)
    # assets today

    aprime = get_aprime(asset_policy, state_index, model_params)
    # assets tomorrow Noption * statesize as assets are contingent on car choice

    Aprime = sum(aprime .* L, dims = 1)[1]

    c = get_consumption(Pces, W, τ_rev, R, asset_policy, model_params, state_index)

    # aggregate asset demand 

    A = sum(a .* L, dims = 1)[1]
    # aggregate asset positoin entering the period

    NetA = -(R * A) + Aprime

    C = sum( c .* L, dims = 1)[1]
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

    production = price_production * sum( TFP .* ef_units .* labor_supply .* L, dims = 1)[1] 

    income = sum( wz.* L, dims = 1)[1] 
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

function make_dataset(output, output_int)

    backfilllength = 10

    bigI = []
    bigT = []
    bigC = []
    bigN = []
    bigIHS_car = []
    bigLFP = []
    country_index = []

    backfill= Array{Float64}(undef, backfilllength )

    for cnt = 1:trade_params().Ncntry
        
        ###############################################
        # Investment 
        fill!(backfill, output_int[cnt].I[1])

        I = [output[cnt][xxx].I[1] for xxx in 1:T]

        prepend!(I, backfill)
    
        time = -9:1:(length(I)-10)
    
        append!(bigI, I)
        append!(bigT,time)
        append!(country_index, Int.(cnt.*ones(length(time))))
    
        ###############################################
        # Counsumption
        fill!(backfill, output_int[cnt].C[1])

        C = [output[cnt][xxx].C[1] for xxx in 1:T]

        prepend!(C, backfill)
        append!(bigC, C)
    
        ###############################################
        # Labor Supply
        fill!(backfill, output_int[cnt].N[1])

        N = [output[cnt][xxx].N[1] for xxx in 1:T]

        prepend!(N, backfill)
        append!(bigN, N)

        ###############################################
        # CARS
        fill!(backfill, output_int[cnt].IHS_car[1])

        IHS_car = [output[cnt][xxx].IHS_car[1] for xxx in 1:T]

        prepend!(IHS_car, backfill)
        append!(bigIHS_car, IHS_car)
        ###############################################
        # LFP
        fill!(backfill, output_int[cnt].LFP[1])

        LFP = [output[cnt][xxx].LFP[1] for xxx in 1:T]

        prepend!(LFP, backfill)
        append!(bigLFP, LFP)
    end

    df = DataFrame(country_index = country_index, 
               consumption = bigC,
               labor_supply = bigN,
                durable_c = bigI,
                time = bigT,
                cars = bigIHS_car,
                LFP = bigLFP
               );

    return df

end