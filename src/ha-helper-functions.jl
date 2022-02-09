struct NIPA{T}
    C::Array{T}
    I::Array{T}
    income::Array{T}
    N::Array{T}
    Aprime::Array{T}
    IHS_car::Array{T}
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

function get_consumption(Pces, W, τ_rev, R, model_params, asset_policy, car_policy, work_policy, state_index)
 
    @unpack Na, Nshocks, Doptions, Woptions, mc, dgrid, agrid,
    pnew, pused, statesize = model_params

    c = Array{eltype(R)}(undef, statesize)
    fill!(c, 0.0)
    
    invest = similar(c)
    fill!(invest,0.0)

    new_car = similar(c)
    fill!(new_car,0.0)

    @inbounds for (foo, xxx) in enumerate(state_index)

        for car = 1:Doptions

            ihold = investment(car, dgrid[xxx[3]], dgrid[end], pnew, pused)

            if car == 1

                new_car_hold = 1.0

            elseif car == 2
                new_car_hold = 0.0

            end

            for wrk = 1:Woptions
            # work through different work options

                wz = labor_income(exp.(mc.state_values[xxx[2]]), W, wrk)
            # will return labor income depending upon how much working.

                share_type = work_policy[xxx[1],xxx[2], xxx[3], car, wrk] * car_policy[xxx[1],xxx[2], xxx[3], car] 
                # share of worker type. work*car
                
                invest[foo] += share_type * ihold

                new_car[foo] += share_type * new_car_hold
                       
                chold = sum( asset_policy[xxx[1], : ,xxx[2], xxx[3], car, wrk] .* consumption(Pces, τ_rev, R*agrid[xxx[1]], agrid, 
                            wz, ihold )) # given different aprim, all the consumption 

                c[foo] += share_type * chold

            end
        end

    end

    return c, invest, new_car

end

##############################################################################

function get_laborsupply(model_params, W, distribution, car_policy, work_policy, state_index)

    @unpack Na, Ncars, Nshocks, Woptions, Doptions, agrid,
    statesize = model_params
       
    wz, ef_units = get_laborincome(W, model_params, car_policy, work_policy, state_index)

    labor_supply = Array{eltype(work_policy)}(undef, statesize)

    fill!(labor_supply, 0.0) #need to fill given += operator below

    for (foo, xxx) in enumerate(state_index)

        for car = 1:Doptions

            labor_supply[foo] += car_policy[xxx[1],xxx[2], xxx[3], car] * work_policy[xxx[1],xxx[2], xxx[3], car, 1]

        end

    end

    LFP = sum(labor_supply .* distribution, dims = 1)[1]

    N = sum( ef_units .* labor_supply .* distribution, dims = 1)[1]

    return labor_supply, N, LFP
    # matrix that is of size Na*Nshocks*Ncars

end

##############################################################################

function get_distribution(state_index, stationary_distribution)
    # returns either the joint (default) distribution, asset if key word is passed,

    # assets is first index,
    # shock is second index,
    # car is thrid index 

    asset_state = [xxx[1] for xxx in state_index]
        # asset should be second index

    distribution = [sum(stationary_distribution[asset_state .== xxx]) for xxx in unique(asset_state)]

    return distribution

end

##############################################################################

function get_laborincome(W, model_params, car_policy, work_policy, state_index)
    # grabs the current shock in level terms
    # asset is first index,
    # shock is second index,
    # car is thrid index 
    @unpack statesize, Woptions, Doptions, mc = model_params
    
    wz = Array{eltype(mc.state_values)}(undef, statesize)
    fill!(wz,0.0)

    ef_units = similar(wz)

    for (foo, xxx) in enumerate(state_index)

        ef_units[foo] = exp.(mc.state_values[xxx[2]])

        for car = 1:Doptions

            for wrk = 1:Woptions
            # work through different work options

                wz[foo] += car_policy[xxx[1],xxx[2], xxx[3], car] * 
                work_policy[xxx[1],xxx[2], xxx[3], car, wrk] * labor_income(exp.(mc.state_values[xxx[2]]), W, wrk)
            # will return labor income depending upon how much working.

            end

        end

    end

    return wz, ef_units

end


##############################################################################

function get_astate(model_params, state_index)
    # grabs the current asset level
    # shock is first index,
    # assets is second index,
    # car is thrid index 
    @unpack statesize, agrid = model_params
    
    asset_state = Array{eltype(agrid)}(undef, statesize, 1)

    for (foo, xxx) in enumerate(state_index)
        
        asset_state[foo, 1] = agrid[xxx[1]]

    end

    return asset_state

end

##############################################################################

function get_expenditure(Pces, W, τ_rev, R, model_params, asset_policy, 
    car_policy, work_policy, state_index, distribution)

    c, invest, new_car = get_consumption(Pces, W, τ_rev, R, model_params, asset_policy, car_policy, work_policy, state_index)

    # aggregate asset demand 

    C = sum( c .* distribution, dims = 1)[1]
    #aggregate nondurable consumption
    # then multiply the prob of being in that state and sum to aggregate 

    I = sum( invest .* distribution, dims = 1)[1]

    a = get_astate(model_params, state_index)
    # assets today

    aprime = get_aprime(model_params, asset_policy, car_policy, work_policy, state_index)
    # assets tomorrow Noption * statesize as assets are contingent on car choice

    Aprime = sum(aprime .* distribution, dims = 1)[1]

    A = sum(a .* distribution, dims = 1)[1]
    # aggregate asset positoin entering the period

    NetA = -(R * A) + Aprime

    return C + I + NetA
    # don't quite undsertand this...path only works/makes sense if this way
    # neet to factor in how these resources are evolving. 

end


##############################################################################

function aggregate(Pces, W, τ_rev, R, TFP, model_params, asset_policy, 
    car_policy, work_policy, state_index, distribution; display = false)

    wz, ef_units = get_laborincome(W, model_params, car_policy, work_policy, state_index)
    #labor income shocks

    a = get_astate(model_params, state_index)
    # assets today

    aprime = get_aprime(model_params, asset_policy, car_policy, work_policy, state_index)
    # assets tomorrow Noption * statesize as assets are contingent on car choice

    Aprime = sum(aprime .* distribution, dims = 1)[1]

    c, invest, new_car = get_consumption(Pces, W, τ_rev, R, model_params, asset_policy, car_policy, work_policy, state_index)

    # aggregate asset demand 

    A = sum(a .* distribution, dims = 1)[1]
    # aggregate asset positoin entering the period

    NetA = -(R * A) + Aprime

    C = sum( c .* distribution, dims = 1)[1]
    #aggregate nondurable consumption
    # then multiply the prob of being in that state and sum to aggregate 

    
    IHS_car = sum( new_car .* distribution, dims = 1)[1]
    #aggregate nondurable consumption
    # then multiply the prob of being in that state and sum to aggregate 

    I = sum( invest .* distribution, dims = 1)[1]
    #aggregate durable expenditure
    # same for durable

    expenditure = Pces *( C + I + NetA - τ_rev)
    # seems neccessary to get everything in same units 
    # when price levels are different

    labor_supply, N, LFP = get_laborsupply(model_params, W, distribution, car_policy, work_policy, state_index)

    price_production = W / TFP # comes of profit max condition that p = w / A

    production = price_production*sum( TFP .* ef_units .* labor_supply .* distribution, dims = 1)[1] 
    # tariff revenue shows up here like production of G, same with income...
    # needs to be here to line up with expenditure.

    income = sum( wz.* distribution, dims = 1)[1] 
    # aggregate labor income 
    # take shock and multiply by prob of being in that state, sum to aggregate

    output_stats = NIPA(
        [C], [I], [income], [N], [Aprime], [IHS_car], [LFP]
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
        println("Aggregate Durable Expenditure")
        println(round(I, digits = digits))
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