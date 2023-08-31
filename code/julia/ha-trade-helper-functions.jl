struct NIPA
    PC::Float64 #consumption
    M::Float64 #imports
    X::Float64 #exports
    income::Float64
    production::Float64
    N::Float64
    Aprime::Float64
    NetA::Float64
    G::Float64
end

struct trade
    bilateral_imports::Array{Float64} # asset_policy
    bilateral_πprob::Array{Float64} # choice probabilities
    tariff_revenue::Float64
end

struct hhXsection
    wz::Array{Float64} #consumption
    a::Array{Float64} #imports
    income::Array{Float64} #exports
    pc::Array{Float64}
    homeshare::Array{Float64}
    mpc_avg::Array{Float64}
    welfare::Array{Float64}
    θavg::Array{Float64}
end

struct micromoments
    poor_πii::Float64 #consumption
    rich_πii::Float64 #imports
    middle_πii::Float64 #exports
    poor_θ::Float64
    rich_θ::Float64
    middle_θ::Float64
end
##############################################################################

function cal_make_stats(Xsec ; prctile = [50.0, 50.0])
    # takes df that is the form of hhXsection

    poor = Xsec.income .< percentile(Xsec.income, prctile[1])
    rich = Xsec.income .> percentile(Xsec.income, prctile[2])

    middle = (Xsec.income .> percentile(Xsec.income, 55.0)) .== (Xsec.income .< percentile(Xsec.income, 45.0))

    poor_πii = median(Xsec.homeshare[poor])
    rich_πii = median(Xsec.homeshare[rich])
    middle_πii = median(Xsec.homeshare[middle])

    poor_θ = median(Xsec.θavg[poor])
    rich_θ = median(Xsec.θavg[rich])
    middle_θ = median(Xsec.θavg[middle])

    return micromoments(poor_πii, rich_πii, middle_πii, poor_θ, rich_θ, middle_θ)

end

function make_stats(df; prctile = [50, 50])
    # takes df that is the form of hhXsection

    poor = df.income .< percentile(df.income, prctile[1])
    rich = df.income .> percentile(df.income, prctile[2])

    middle = (df.income .> percentile(df.income, 55)) .== (df.income .< percentile(df.income, 45))

    poor_πii = median(df[poor,:].homeshare)
    rich_πii = median(df[rich,:].homeshare)
    middle_πii = median(df[middle,:].homeshare)

    poor_θ = median(df[poor,:].θ)
    rich_θ = median(df[rich,:].θ)
    middle_θ = median(df[middle,:].θ)

    poor_mpc = median(df[poor,:].mpc)
    rich_mpc = median(df[rich,:].mpc)
    middle_mpc= median(df[middle,:].mpc)

    poor_∂W = median(df[poor,:].∂W)
    rich_∂W = median(df[rich,:].∂W)
    middle_∂W= median(df[middle,:].∂W)

    return (rich_πii, rich_θ, rich_∂W, rich_mpc), (poor_πii, poor_θ, poor_∂W, poor_mpc), 
        (middle_πii, middle_θ, middle_∂W ,middle_mpc)

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

function get_consumption(p, cons_policy, πprob, state_index, model_params)
    # grabs the choice of assets given states

    @unpack Na, Nshocks, Ncntry = model_params
    
    pconsumption = Array{eltype(cons_policy)}(undef, Na*Nshocks)
    fill!(pconsumption, 0.0) #need to fill given += operator below

    for (foo, xxx) in enumerate(state_index)

        for cntry = 1:Ncntry

            pconsumption[foo] +=  p[cntry] * cons_policy[xxx[1], xxx[2], cntry] * πprob[xxx[1], xxx[2], cntry]

        end
        
    end

    return pconsumption

end

##############################################################################

function get_trade(p, cons_policy, πprob, state_index, model_params)
    # grabs the choice of assets given states

    @unpack Na, Nshocks, Ncntry = model_params
    
    pc_by_state = Array{eltype(cons_policy)}(undef, Na*Nshocks, Ncntry)
    pcπ_by_state = Array{eltype(cons_policy)}(undef, Na*Nshocks, Ncntry)
    πprob_by_state = Array{eltype(cons_policy)}(undef, Na*Nshocks, Ncntry)

    fill!(pc_by_state , 0.0) #need to fill given += operator below
    fill!(pcπ_by_state , 0.0) #need to fill given += operator below
    fill!(πprob_by_state, 0.0) #need to fill given += operator below

    for (foo, xxx) in enumerate(state_index)

        for cntry = 1:Ncntry

            #pc = ( -asset_policy[xxx[1], xxx[2], cntry] + R*agrid[xxx[1]] + W*shocks[xxx[2]] )
            #should watch this line

            pc = p[cntry] * cons_policy[xxx[1], xxx[2], cntry]

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

function aggregate(R, W, p, τ, tariff, country, household, distribution, hh_params; display = false)

    # organization...

    @unpack state_index, λ = distribution
    @unpack πprob, asset_policy, cons_policy = household
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

    pc = get_consumption(p, cons_policy, πprob, state_index, hh_params)

    pcπ_by_state, π_by_state = get_trade(p, cons_policy, πprob, state_index, hh_params)[1:2]

    Aprime = L .* dot(aprime, λ) 
    # asset holdings next period (17)

    # aggregate asset demand 

    A = L * dot(a , λ)
    # aggregate asset positoin entering the period

    NetA = -(R * A) + Aprime

    PC = L * dot(pc , λ)
    #aggregate consumption

    G = L * τ

    ##############
    # the Trade Stuff

    imports = sum(pcπ_by_state[:, 1:end .!= country], dims = 2)
    # exclude home country imports

    M = L * dot(imports , λ)
    # aggregate imports by value

    bilateral_imports = sum( L * pcπ_by_state.* λ, dims = 1)

    tariff_revenue  =  sum(( tariff[country, :] ./ ( 1.0 .+ tariff[country, :]) ).* bilateral_imports' )

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
        println("Government Spending")
        println(round(G, digits = digits))
        println("Tariff Revenue")
        println(tariff_revenue)

    end

    return NIPA(PC, M, X, income, production, N, Aprime, NetA, G), trade(bilateral_imports, bilateral_πprob, tariff_revenue)

end

##############################################################################
##############################################################################

function make_Xsection(R, W, p, household, distribution, θ, mpc, ∂W, home_cntry, model_params; Nsims = 100)
    
    @unpack mc, agrid, Ncntry = model_params
    @unpack state_index = distribution
    @unpack cons_policy, πprob = household
    @unpack θπ, θπii, θc, θcii = θ

    Qmc = MarkovChain(distribution.Q) 
    # converts markov chain to a markov chain interperted by 
    # quant econ
    Random.seed!(03281978) # this sets the seed
    X = state_index[simulate(Qmc, Nsims)]
    # this returns the states

    a = Array{eltype(W)}(undef, Nsims)

    homeshare = similar(a)

    θavg = similar(a)

    income = similar(a)

    mpc_avg = similar(a)

    welfare = similar(a)

    θx = Array{eltype(W)}(undef, Ncntry-1)

    θweight = similar(θx)

    weight = similar(θx)

    mpc_weight = Array{eltype(W)}(undef, Ncntry)

    pc = similar(a)

    fill!(pc, 0.0)

    ef_units = exp.(mc.state_values)

    # cons_πprob = cons_policy .* πprob

    @inbounds for (foo, xxx) in enumerate(X)

        a[foo] = agrid[xxx[1]]

        income[foo] = labor_income(ef_units[xxx[2]], W)
        # + R*a[foo]

        cntry_count = 0

        # fill!(weight, 0.0)
        # fill!(θweight, 0.0)
        # fill!(mpc_weight, 0.0)

        for cntry = 1:Ncntry

            cntry_pc = p[cntry] * cons_policy[xxx[1], xxx[2], cntry] * πprob[xxx[1], xxx[2], cntry]          
            # 

            pc[foo] +=  cntry_pc

            mpc_weight[cntry] = mpc[xxx[1], xxx[2], cntry]*cntry_pc

            if cntry != home_cntry

                cntry_count = cntry_count + one(cntry_count)

                weight[cntry_count] = cntry_pc

                θweight[cntry_count] = (1.0 + ( θπ[xxx[1], xxx[2], cntry] - θπii[xxx[1], xxx[2],cntry]) + 
                ( θc[xxx[1], xxx[2], cntry] - θcii[xxx[1], xxx[2], cntry] ) ) * weight[cntry_count]

            end
            
        end

        welfare[foo] = ∂W[xxx[1], xxx[2]] 

        homeshare[foo] = ( p[home_cntry] * cons_policy[xxx[1], xxx[2], home_cntry] 
                    * πprob[xxx[1], xxx[2], home_cntry] ) / pc[foo]

        θavg[foo] = sum(θweight) / sum(weight)

        mpc_avg[foo] = sum(mpc_weight ) / pc[foo]

    end

    return hhXsection(income, a, income, pc, homeshare, mpc_avg, welfare, θavg)

end




##############################################################################

function bilateral_consumption(R, W, hh, country, model_params)
    # grabs the choice of assets given states

    @unpack Na, Nshocks, Ncntry, mc, agrid = model_params
    
    pconsumption = Array{Float64}(undef, Na*Nshocks, Ncntry)

    avg_expenditure = Array{Float64}(undef, Na, Nshocks)
    fill!(avg_expenditure, 0.0)

    shocks = exp.(mc.state_values)

    state_index = Array{Tuple{eltype(Int64), eltype(Int64)}}(undef, Na*Nshocks, 1)
    
    make_state_index!(state_index, model_params)

    for (foo, xxx) in enumerate(state_index)

        for cntry = 1:Ncntry

            pconsumption[foo, cntry] =  ( -hh[country].asset_policy[xxx[1], xxx[2], cntry] 
                    + R[country]*agrid[xxx[1]] + W[country]*shocks[xxx[2]] ) * hh[country].πprob[xxx[1], xxx[2], cntry]
                    
            avg_expenditure[xxx[1], xxx[2]] += hh[country].πprob[xxx[1], xxx[2], cntry] .* ( -hh[country].asset_policy[xxx[1], xxx[2], cntry] 
            + R[country]*agrid[xxx[1]] + W[country]*shocks[xxx[2]] )
            
            # hh[country].cons_policy[xxx[1], xxx[2], cntry]

        end
        
    end

    return pconsumption, avg_expenditure
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
    #todo need to add other cases and incorperate

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
