function eq_variation(R, w, p, Δ_hh, state_index, model_params)
    # these should be the old prices
    # Δ_hh is the value function at the new prices
    # this works through everything state by state


    @unpack σϵ, Na, Nshocks, Ncntry = model_params

    τeqv = Array{Float64}(undef,Na,Nshocks)

    Δ_v = log_sum_column(Δ_hh.Tv,  σϵ)
    
    xguess = [0.0]

    n = length(xguess)
    diag_adjust = n - 1

    for (foo, xxx) in enumerate(state_index)
        # work through all the states

        # ΔW[xxx[1], xxx[2]] = eq_variation(τ , xxx[1], xxx[2], R, w, p, Δ_v[xxx[1], xxx[2]], model_params)

        f(x) = eq_variation(x, xxx[1], xxx[2], R, w, p, Δ_v[xxx[1], xxx[2]], model_params)

        function f!(fvec, x)

            fvec .= f(x)
        
        end

        # find the transfer to make a guy indifferent
        sol = fsolve(f!, xguess, show_trace = true, method = :hybr;
            ml=diag_adjust, mu=diag_adjust,
            diag=ones(n),
            mode= 1,
            tol=1e-10,)


        τeqv[xxx[1], xxx[2]] = sol.x[1]

    end
    
    return τeqv

end

#########################################################################################

function eq_variation(xxx, astate, shockstate, R, w, p, Δ_v, model_params)
    # core function that computes difference between 
    # value fun at old prices + transfer (xxx) and new value fun

    eqv_hh = solve_household_problem(R, w, p, xxx[1], model_params)

    v = log_sum_v(eqv_hh.Tv[astate, shockstate,:], model_params.σϵ, model_params.Ncntry)

    return Δ_v - v
    #Δ_v is new value fun

end

##############################################################################
##############################################################################

function lucas_eq_variation(xxx, astate, shockstate, hh, Δhh, model_params)
    # core function that computes difference between 
    # value fun at old prices + transfer (xxx) and new value fun

    @unpack ψ , σϵ , Ncntry = model_params

    Tv = value_function_fixedpolicy(hh, xxx[1], model_params)    

    v = log_sum_v(ψ[astate, shockstate, :] .+ Tv[astate, shockstate, : ], σϵ, Ncntry)

    Δ_v = log_sum_v(ψ[astate, shockstate, :] .+ Δhh.Tv[astate, shockstate, : ], σϵ, Ncntry)

    return Δ_v - v
    #Δ_v is new value fun

end

##############################################################################

function lucas_eq_variation(hh, Δhh, state_index, model_params)
    # hh should be the old prices
    # Δ_hh is the value function at the new prices
    # this works through everything state by state


    @unpack σϵ, Na, Nshocks = model_params

    λeqv = Array{Float64}(undef, Na, Nshocks)

    xguess = [1.0]

    n = length(xguess)
    diag_adjust = n - 1

    for (foo, xxx) in enumerate(state_index)
        # work through all the states

        # ΔW[xxx[1], xxx[2]] = eq_variation(τ , xxx[1], xxx[2], R, w, p, Δ_v[xxx[1], xxx[2]], model_params)

        f(x) = lucas_eq_variation(x, xxx[1], xxx[2], hh, Δhh, model_params)

        function f!(fvec, x)

            fvec .= f(x)
        
        end

        # find the transfer to make a guy indifferent
        sol = fsolve(f!, xguess, show_trace = false, method = :hybr;
            ml=diag_adjust, mu=diag_adjust,
            diag=ones(n),
            mode= 1,
            tol=1e-10,)


        λeqv[xxx[1], xxx[2]] = sol.x[1]

    end
    
    return λeqv

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

function welfare_by_state(hh, Δ_hh, dist, Δ_dist, country, model_params)

    @unpack σϵ, ψ = model_params

    v = alt_log_sum(hh[country].πprob, hh[country].Tv, σϵ, ψ)    

    Δ_v = alt_log_sum(Δ_hh[country].πprob, Δ_hh[country].Tv, σϵ, ψ)

    SW = sum(vec(v).*dist[country].λ)

    Δ_SW = sum(vec(Δ_v).*Δ_dist[country].λ)

    ∂SW = ( Δ_SW - SW) / SW

    ∂W = ( Δ_v .- v) ./ SW

    return ∂W, ∂SW, v , Δ_v

end

##############################################################################

function value_function_fixedpolicy(household, λ, model_params; tol = 10^-10, Niter = 500)
    # If you give me policy functions, 
    # I construct the value function
    
    @unpack Na, Nshocks, Ncntry, β, σϵ, ψ = model_params
    @unpack cons_policy, asset_policy, πprob = household

    # this is the guess... always start at borrowing cosntraint
    v = -ones(Na, Nshocks, Ncntry)/(1-β)
    Tv = similar(v)

    for iter in 1:Niter
        
        make_Tv!(Tv, v, cons_policy, asset_policy, πprob, λ, ψ, model_params)

        err = vec_max(Tv, v)

        #println(iter)

        if err < tol

            #println(iter)

            break
        end

        copy!(v,Tv)

        if iter == Niter

          println("value function may not have converged")
          println("check the situation")
          
        end

    end

    return Tv
    
end