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

function get_asset_shock_state(state_index, model_params)

    @unpack Na, Nshocks, agrid = model_params

    asset_state = Array{Float64}(undef, Na*Nshocks)

    shock_state = Array{Float64}(undef, Na*Nshocks)

    for (foo, xxx) in enumerate(state_index)
        
        asset_state[foo] = agrid[xxx[1]]

        shock_state[foo] = xxx[2]

    end

    return asset_state, shock_state

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

function aggregate(R, W, p, country, household, distribution, hh_params; display = false)

    # organization...

    @unpack state_index, λ = distribution
    @unpack πprob, asset_policy = household
    @unpack TFP, L = hh_params

    #####
    wz, ef_units = get_laborincome(W, state_index, hh_params)

    N = L* dot(ef_units, λ) #number of effeciency units (16)

    #####
    # Get stuff from hh side

    a = get_astate(state_index, hh_params)
    # assets today

    aprime = get_aprime(asset_policy, πprob, state_index, hh_params)
    # # assets tomorrow Noption * statesize as assets are contingent on car choice

    pc = get_consumption(R, W, asset_policy, πprob, state_index, hh_params)

    pcπ_by_state, π_by_state = get_trade(R, W, asset_policy, πprob, state_index, hh_params)[1:2]

    Aprime = L .* dot(aprime, λ) 
    # asset holdings next period (17)

    # aggregate asset demand 

    A = L * dot(a , λ)
    # aggregate asset positoin entering the period

    NetA = -(R * A) + Aprime

    PC = L * dot(pc , λ)
    #aggregate consumption

    ##############
    # the Trade Stuff

    imports = sum(pcπ_by_state[:, 1:end .!= country], dims = 2)
    # exclude home country imports

    M = L * dot(imports , λ)
    # aggregate imports by value

    bilateral_imports = sum( L * pcπ_by_state.* λ, dims = 1)

    bilateral_πprob = sum( π_by_state.* λ, dims = 1)

    production = p[country] * TFP* N 

    home_consumption = sum(pcπ_by_state[:, 1:end .== country], dims = 2)

    X = production - L * dot(home_consumption , λ)

    #############
    # then compute GDP like measure

    income = L * dot(wz, λ)

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
function social_welfare(hh, Δ_hh, dist, Δ_dist, country, σϵ)

    W = sum(vec(log_sum_column(hh[country].Tv,  σϵ)).*dist[country].λ)

    Δ_W = sum(vec(log_sum_column(Δ_hh[country].Tv,  σϵ)).*Δ_dist[country].λ)

    ∂W = 100.0 .*( W - Δ_W) ./ W  

    lost = vec(log_sum_column(hh[country].Tv,  σϵ)) .> vec(log_sum_column(Δ_hh[country].Tv,  σϵ))

    share_lost = sum(dist[country].λ[lost])
    
    return ∂W, share_lost

end

##############################################################################

function welfare_by_state(hh, Δ_hh, country, σϵ)

    v = log_sum_column(hh[country].Tv,  σϵ)

    Δ_v = log_sum_column(Δ_hh[country].Tv,  σϵ)

    ∂logW = 100.0 .*( v - Δ_v) ./ v  
    
    ∂W = ( v - Δ_v)

    return ∂W, ∂logW

end

##############################################################################

function bilateral_consumption(R, W, hh, country, model_params)
    # grabs the choice of assets given states

    @unpack Na, Nshocks, Ncntry, mc, agrid = model_params
    
    pconsumption = Array{Float64}(undef, Na*Nshocks, Ncntry)

    shocks = exp.(mc.state_values)

    state_index = Array{Tuple{eltype(Int64), eltype(Int64)}}(undef, Na*Nshocks, 1)
    
    make_state_index!(state_index, model_params)

    for (foo, xxx) in enumerate(state_index)

        for cntry = 1:Ncntry

            pconsumption[foo, cntry] =  ( -hh[country].asset_policy[xxx[1], xxx[2], cntry] 
                    + R[country]*agrid[xxx[1]] + W[country]*shocks[xxx[2]] ) * hh[country].πprob[xxx[1], xxx[2], cntry]

        end
        
    end

    return pconsumption

end

##############################################################################

function πii_elasticity(hh, Δ_hh, d, Δ_d, country, model_params)

    @unpack Na, Nshocks, Ncntry, mc, agrid = model_params

    ∂log_πii = Array{Float64}(undef, Na*Nshocks, 1)

    ∂log_πij = Array{Float64}(undef, Na*Nshocks, Ncntry)

    πij = Array{Float64}(undef, Na*Nshocks, Ncntry)

    state_index = Array{Tuple{eltype(Int64), eltype(Int64)}}(undef, Na*Nshocks, 1)
    
    make_state_index!(state_index, model_params)

    τ_change = median(log.(Δ_d[:,country]) .- log.(d[:,country]))

    for (foo, xxx) in enumerate(state_index)
        
        ∂log_πii[foo] = (log.(Δ_hh[country].πprob[xxx[1],xxx[2], country] ) 
            .- log.(hh[country].πprob[xxx[1],xxx[2], country] )) ./ τ_change
        #./ (log.(Δ_d[:,country]) .- log.(d[:,country]))

        for cntry = 1:Ncntry

            ∂log_πij[foo, cntry] =  (log.(Δ_hh[country].πprob[xxx[1],xxx[2], cntry] ) 
            .- log.(hh[country].πprob[xxx[1],xxx[2], cntry] )) ./ ( τ_change )

            πij[foo, cntry] = Δ_hh[country].πprob[xxx[1],xxx[2], cntry] / Δ_hh[country].πprob[xxx[1],xxx[2], country]
        
        end

    end
    
    return ∂log_πii, ∂log_πij, πij

end

##############################################################################

function θ_by_state(pconsumption, Δ_pconsumption, d, Δ_d, country, model_params)

    @unpack Na, Nshocks, Ncntry, mc, agrid = model_params

    θ_micro = Array{Float64}(undef, Na*Nshocks, Ncntry)

    state_index = Array{Tuple{eltype(Int64), eltype(Int64)}}(undef, Na*Nshocks, 1)
    
    make_state_index!(state_index, model_params)

    τ_change = (log.(Δ_d[:,country]) .- log.(d[:,country]))

    for (foo, xxx) in enumerate(state_index)

        trade =  pconsumption[foo, :] ./ pconsumption[foo, country]

        Δ_trade =  Δ_pconsumption[foo, :] ./ Δ_pconsumption[foo, country]
        
        θ_micro[foo, :] = (log.(Δ_trade) .- log.(trade)) ./ τ_change
        #./ (log.(Δ_d[:,country]) .- log.(d[:,country]))

    end
    
    return θ_micro

end

##############################################################################

function make_hh_dataframe(dist, hh, country, R, W, hh_params)

    @unpack Na, Nshocks, mc, agrid = hh_params
    @unpack asset_policy,  πprob = hh[country]
            
    income = Array{Float64}(undef, Na*Nshocks)
    weights = Array{Float64}(undef, Na*Nshocks)
    homeshare = Array{Float64}(undef, Na*Nshocks)
    pconsumption = Array{Float64}(undef, Na*Nshocks)
    home_π = Array{Float64}(undef, Na*Nshocks)

    fill!(income, 0.0) #need to fill given += operator below
    fill!(weights, 0.0) #need to fill given += operator below
    fill!(pconsumption, 0.0) #need to fill given += operator below

    shocks = exp.(mc.state_values)

    for (foo, xxx) in enumerate(dist[country].state_index)

        income[foo] = (1 - R[country])*agrid[xxx[1]] + W[country]*shocks[xxx[2]] 

        weights[foo] = dist[country].λ[foo]

        for cntry = 1:Ncntry

            pconsumption[foo] +=  ( -asset_policy[xxx[1], xxx[2], cntry] 
                    + R[country]*agrid[xxx[1]] + W[country]*shocks[xxx[2]] ) * πprob[xxx[1], xxx[2], cntry]

        end

        home_π[foo] = hh[country].πprob[xxx[1], xxx[2], country]

        homeshare[foo] = (( - asset_policy[xxx[1], xxx[2], country] 
                    + R[country]*agrid[xxx[1]] + W[country]*shocks[xxx[2]] ) 
                    * πprob[xxx[1], xxx[2], country]) ./ pconsumption[foo]

    end
    
    df = DataFrame(income = income, 
               weights = weights,
               homeshare = homeshare,
               home_choice = home_π,
               consumption = pconsumption,
               );
    
    return df
    
end 

##############################################################################

function make_welfare_dataframe(∂W, ∂logW, hh_params)

    @unpack Na, Nshocks, agrid = hh_params
    
    welfare = Array{eltype(∂W)}(undef, Na*Nshocks)
    welfare_level = Array{eltype(∂W)}(undef, Na*Nshocks)
    shock = Array{eltype(∂W)}(undef, Na*Nshocks)
    asset = Array{eltype(∂W)}(undef, Na*Nshocks)
    
    state_index = Array{Tuple{eltype(Int64), eltype(Int64)}}(undef, Na*Nshocks, 1)
    
    make_state_index!(state_index, hh_params)
    
    for (foo, xxx) in enumerate(state_index)
        
        shock[foo] = xxx[2]
        
        asset[foo] = agrid[xxx[1]]

        welfare[foo] = ∂logW[foo]

        welfare_level[foo] = ∂W[foo]
    end
    
    df = DataFrame(asset = asset, 
               shock = shock,
               welfare = welfare,
               welfare_level = welfare_level,
               );
    
    return df
    
end  


function unpack_xvec(xvec, Ncntry)

    if length(xvec) ≈ ( (Ncntry - 1) + Ncntry + 1)
        # this is the financial globalization case
        
        W = [xvec[1:(Ncntry - one(Ncntry))]; 1.0 ]

        τ = xvec[Ncntry:end - 1]

        R = ones(Ncntry)*xvec[end]

    else
        W = NaN

        τ = NaN

        R = NaN

    end

    return W, τ, R

end
