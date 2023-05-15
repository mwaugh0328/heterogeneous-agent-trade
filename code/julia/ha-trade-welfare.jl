function eq_variation(R, w, p, Δ_hh, state_index, model_params)
    # these should be the old prices
    # Δ_hh is the value function at the new prices
    # this works through everything state by state


    @unpack σϵ, Na, Nshocks, Ncntry = model_params

    τeqv = Array{Float64}(undef,Na,Nshocks)

    Δ_v = log_sum_column(Δ_hh.Tv,  σϵ)
    
    xguess = [0.0]

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

    eqv_hh = solve_household_problem(R, w, p, xxx, model_params)

    v = log_sum_v(eqv_hh.Tv[astate, shockstate,:], model_params.σϵ, model_params.Ncntry)

    return Δ_v - v
    #Δ_v is new value fun

end