function transition_path(τ_path, W, τ_rev, R, dist, hh, model_params)
    # main version to be used in solving for eq.
    
    @unpack Ncntry = trade_params()
    
    T == size(τ_path)[1]
    # \tau_path should be an array of arrays...
    # size of it should be 

    @assert T == Int(size(W)[1] / Ncntry)
    
    # then this just organizes things so stuff will
    # be of size Ncntry, T
    W = reshape(W, Ncntry, T)
    τ_rev = reshape(τ_rev, Ncntry, T)
    R = reshape(R, Ncntry, T)

    Pces = similar(W)

    AD = similar(W)
    N = similar(W)

    net_demand = Array{eltype(W)}(undef, 2*Ncntry, T)
    # The output vectov, net demand is 2*Ncountry bc we are 
    # checking if output and tariff revenue is consistent.

    for fwdate = 1:T

        Pces[:, fwdate] = goods_prices(W[:, fwdate], trade_params(τ = τ_path[fwdate]))[2]

    end

    # then do the forward/backward thing for each country
    for nct = 1:Ncntry 

        AD[nct, :], N[nct, :] = foward_backward(Pces[nct,:], W[nct,:], 
                                    τ_rev[nct,:], R[nct,:], dist[nct], hh[nct], model_params) 
        #Tv, L, Pces_int all need to be country specific initial conditions.
        # takes in path of prices, wages, tariff revenue,
    end

    # Then given how each country behaves, we can bring it together through the 
    # trade eq....
    for fwdate = 1:T

        net_demand[:, fwdate] = trade_equilibrium(W[:, fwdate], AD[:, fwdate], N[:, fwdate], τ_rev[:, fwdate], 
                                                    trade_params(τ = τ_path[fwdate]) )

    end

    return net_demand'[:]

end

#####################################################################################################
#####################################################################################################


function foward_backward(Pces, W, τ_rev, R, dist, hh_end, model_params)
    # Main version, if extra arg is passed then it generates output

    @unpack Na, Nshocks, Woptions, mc, β, σa, σw = model_params

    T = length(W)

    hh = Array{household{eltype(W)}}(undef, T, 1)
    output_stats = Array{NIPA{eltype(W)}}(undef, T, 1)

    # setting this up so it's an Array of Array's a bit cleaner, might
    # be more memory effecient?

    AD = Array{eltype(W)}(undef, T)
    N = similar(AD)

    u = Array{eltype(W),4}(undef, Na, Na, Nshocks, Woptions)

    Tv = copy(hh_end.Tv)

    #####################################################################################################
    # This is the backward step: sovle hh problem at T then use bellman operator to work backwards

    @inbounds for bwdate = T:-1:1
        # R is fixed, but it's value depends upon price level last period, that's the value of stuff
        # invested. So you need to take it into account.

        make_utility!(u, Pces[bwdate], W[bwdate], τ_rev[bwdate], R[bwdate], model_params) 

        hh[bwdate] = bellman_operator_policy(Tv, u , mc.p, β, σa, σw);
        
        # use the relationship where
        # Tv_{date} = max [u_{date} + BETV_{date + 1} ]
        # and policy functions for date T are the arg max.

    end

    #####################################################################################################
    # This is the foward step: Start from the inital distribution, L_prime_{1}, then make the Q, given the
    # policy functions at date 1 (from above), gives L_prime_{2} and condtinue forward....
    
    @inbounds for fwdate = 1:T
        # so when date > T as we run it out, just grab stuff from end in policy functions or parameter

        output_stats[xxx] = aggregate(Pces[fwdate], W[fwdate], τ_revenue[fwdate], R[fwdate],
                                     hh[fwdate], dist, A[fwdate], model_params)

        push_foward!(dist, hh[fwdate] , model_params)

    end

    AD = [output[xxx].AD[1] for xxx in 1:T]
    N  = [output[xxx].N[1] for xxx in 1:T]

    return AD, N

end

#####################################################################################################
#####################################################################################################

function push_foward!(dist, household, model_params)
# Pushes the economy forward
# (1) Take policy functions and a distribution L -> aggregates today.
# Policy functions -> Transition probability Q
# Distribution today + Q -> Distribution tomorrow. 

# make transition matrix

make_Q!(dist.Q, dist.state_index, household.asset_policy, household.work_policy, model_params)

# Then push the distribution forward
dist.L = law_of_motion(dist.L , transpose(dist.Q))

end

#####################################################################################################
#####################################################################################################

function collect_end_conditions(W, τ_rev, R, model_params, trade_params)
    # grabs initial conditions for 
    # to compute 

    @unpack Na, Nshocks, Ncars = model_params
    @unpack Ncntry = trade_params

    Pces = goods_prices(W, trade_params)[2]
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    TVend = Array{Array{eltype(W),3}}(undef, Ncntry)

    for ncnty = 1:Ncntry  # do this for each country.
                            # this is the place to use distributed. hh problems can be solved independtly
        TVend[ncnty] = solve_household_problem(Pces[ncnty], W[ncnty], τ_rev[ncnty], R,  model_params)[1]
        
    end

    return TVend

end

##########################################################################
##########################################################################


function collect_intial_conditions(W, τ_revenue, R, model_params, trade_params)
    # grabs initial conditions for 
    # to compute 

    @unpack Na, statesize = model_params
    @unpack Ncntry = trade_params

    Pces = goods_prices(W, trade_params)[2]
    # Given the wage, we need to know the price index to solve the 
    # hh problem. So do this first.

    L = Array{Array{eltype(W)}}(undef, Ncntry)
    
    state_index = Array{Tuple{eltype(Na), 
    eltype(Na), eltype(Na)}}(undef, statesize, 1)

    Q = Array{eltype(W)}(undef, statesize, statesize)

    for ncnty = 1:Ncntry  # do this for each country.

        Q, state_index, L[ncnty] = compute_eq(Pces[ncnty], W[ncnty], τ_revenue[ncnty], R,  model_params)[1:3];

    end

    return Q, state_index, L

end



