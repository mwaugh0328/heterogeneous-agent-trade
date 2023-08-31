#functions used to calibrate the model


function calibrate(xxx, R, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")
    # multiple dispatch version that creates own initial guess

    @unpack Ncntry = cntry_params

    dfguess = DataFrame(guess = xxx);

    CSV.write("current-guess.csv", dfguess)

    ## A bunch of organization here ####################

    @assert length(xxx) == 2*(Ncntry - 1) + 4 + 6  

    # TFP, d = make_country_params(xxx, cntry_params, grvparams, trade_cost_type = trade_cost_type)

    dfWguess = DataFrame(CSV.File("wage-guess.csv"))

    initial_x = log.([dfWguess.guess[1:18]; 0.93])

    return calibrate(xxx, R, initial_x, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
            hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")

end


function calibrate(xxx, R, initial_x, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")
    # this takes primitives (i) solves for an eq. and then (ii) runs gravity regression

    @unpack Ncntry = cntry_params

    dfguess = DataFrame(guess = xxx);

    CSV.write("current-guess.csv", dfguess)

    ## A bunch of organization here ####################

    @assert length(xxx) == 2*(Ncntry - 1) + 4 + 6  

    TFP, d = make_country_params(xxx, cntry_params, grvparams, trade_cost_type = trade_cost_type)

    ##################################################################

    calibrate_cntry_params = country_params(TFP = TFP, d = d, 
                            Ncntry = Ncntry, L = cntry_params.L)


    f(x) = world_equillibrium_FG(exp.(x), R, hh_params, calibrate_cntry_params; tol_vfi = tol_vfi, tol_dis = tol_dis, 
    hh_solution_method = hh_solution_method, stdist_sol_method=stdist_sol_method);

    function f!(fvec, x)
    
        fvec .= f(x)
    
    end
        
    n = length(initial_x)
    diag_adjust = n - 1
    
    sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
          ml=diag_adjust, mu=diag_adjust,
          diag=ones(n),
          mode= 1,
          tol=tol_vfi,
           )
    
    #print(sol)
    
    Wsol = [exp.(sol.x[1:(Ncntry - 1)]); 1.0]
    Wsol = Wsol ./ ( sum(Wsol / Ncntry) )

    
    Rsol = ones(Ncntry)*R

    β = exp(sol.x[end])

    foo_hh_params = household_params(hh_params, β = β)

    πshare= world_equillibrium(Rsol, Wsol, foo_hh_params, calibrate_cntry_params)[5];
    # using the mulitiple dispacth version...since calibration tariffs are zero

    ##################################################################
    # Run gravity regression on model "data"

    trademodel = log.(normalize_by_home_trade(πshare, Ncntry)')

    dfmodel = hcat(DataFrame(trade = vec(drop_diagonal(trademodel, Ncntry))), grvparams.dfcntryfix)

    grvmodel = gravity(dfmodel, trade_cost_type =  trade_cost_type, display = true)

    out_moment_vec = [grvmodel.S[1:end-1] .- grvdata.S[1:end-1] ; 
        grvmodel.θm[1:end-1] .- grvdata.θm[1:end-1] ;
        grvmodel.dist_coef .- grvdata.dist_coef;
        grvmodel.lang_coef .- grvdata.lang_coef]

    ##################################################################

    return out_moment_vec, Wsol, β, πshare
end

# #####################################################################################################

# function calibrate_micro(xxx, micro_moments, R, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
#     hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")
#     # multiple dispatch version that creates own initial guess

#     @unpack Ncntry = cntry_params

#     dfguess = DataFrame(guess = xxx);

#     CSV.write("current-guess.csv", dfguess)

#     ## A bunch of organization here ####################

#     @assert length(xxx) == 2*(Ncntry - 1) + 4 + 6  

#     # TFP, d = make_country_params(xxx, cntry_params, grvparams, trade_cost_type = trade_cost_type)

#     dfWguess = DataFrame(CSV.File("wage-guess.csv"))

#     initial_x = log.([dfWguess.guess[1:18]; 0.93])

#     return calibrate_micro(xxx, micro_moments, R, initial_x, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
#             hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")

# end

# #####################################################################################################

# function calibrate_micro(xxx, micro_moments, R, initial_x, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
#     hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")

#     out_moment_vec, Wsol, β, πshare, hh, dist = calibrate(xxx, R, initial_x, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
#             hh_solution_method = hh_solution_method, stdist_sol_method = stdist_sol_method , trade_cost_type = trade_cost_type)

#     foo_hh_params = household_params(hh_params, β = β)

#     TFP, d = make_country_params(xxx, cntry_params, grvparams, trade_cost_type = trade_cost_type)

#     @unpack σϵ, γ, ψslope = foo_hh_params

#     cntry = 19

#     p = make_p(Wsol, TFP, d[cntry, :], cntry_params.tariff[cntry, :] )

#     ψ = make_ψ(cntry, ψslope.*TFP[cntry].^(1.0 - γ), foo_hh_params)

#     agrid = make_agrid(hh_params, TFP[cntry])

#     cntry_hh_prm = household_params(foo_hh_params, agrid = agrid, 
#     TFP = TFP[cntry], L = L[cntry], σϵ = σϵ*(TFP[cntry]^(1.0 - γ)), ψ = ψ)

#     println("here")

#     θ = make_θ(cntry, Rsol[cntry], Wsol[cntry], p, 0.0, cntry_hh_prm; points = 3, order = 1)
      # line keeps triggering a segfault...no idea why  

#     τeqv = zeros(cntry_hh_prm.Na, cntry_hh_prm.Nshocks)

#     fooX = make_Xsection(Rsol[cntry], Wsol[cntry], p, hh[cntry], dist[cntry],
#          θ, τeqv, τeqv, cntry, cntry_hh_prm; Nsims = 100000)

#     # # println("here")     

#     # # microm = cal_make_stats(fooX)

#     return out_moment_vec, Wsol, β,  πshare, fooX 

# end


#####################################################################################################

function make_country_params(xxx, cntry_params, gravity_params; trade_cost_type = "ek")

    @unpack Ncntry = cntry_params

    TFP = [exp.(xxx[1:(Ncntry - 1)]); 1.0] # S's are normalized -> only have 18 degrees of freedom on Ts
    
    θm = [xxx[Ncntry:((Ncntry - 1)*2)]; -sum(xxx[Ncntry:((Ncntry - 1)*2)])] # same with this, they sum zero

    dist_coef = xxx[((Ncntry - 1)*2 + 1):((Ncntry - 1)*2 + 6)] # six distance bins
    
    lang_coef = xxx[((Ncntry - 1)*2 + 7):end] # the language stuff

    d = zeros(Ncntry,Ncntry)

    make_trade_costs!(trade_costs(dist_coef, lang_coef, θm), d, gravity_params, trade_cost_type = trade_cost_type)

    return TFP, d

end

##########################################################################
##########################################################################

# function calibrate_world_equillibrium(xxx, R, grvdata, grv_params, hh_params, cntry_params; tol_vfi = 1e-10, tol_dis = 1e-10, 
#     hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")

# # This is the function to try and do everything in one step, 
# # still not working...

#     @unpack Ncntry = cntry_params

#     ## A bunch of organization here ####################

#     prices = xxx[1:Ncntry]

#     W = [prices[1:(Ncntry - one(Ncntry))]; 1.0]

#     R = ones(Ncntry)*R

#     foo_hh_params = household_params(hh_params, β = xxx[end])

#     TFP_grav_params = xxx[Ncntry+1:end]

#     @assert length(TFP_grav_params) == 2*(Ncntry - 1) + 4 + 6  

#     TFP, d = make_country_params(TFP_grav_params, cntry_params, grv_params, trade_cost_type = trade_cost_type)

#     ##################################################################

#     calibrate_cntry_params = country_params(TFP = TFP, d = d, Ncntry = Ncntry, L = cntry_params.L)

#     Y, tradeflows, A_demand, Gbudget, πshare = world_equillibrium(R, W, foo_hh_params, calibrate_cntry_params, tol_vfi = tol_vfi, tol_dis = tol_dis, 
#         hh_solution_method = hh_solution_method, stdist_sol_method=stdist_sol_method)[1:5]

#     goods_market = Y .- vec(sum(tradeflows, dims = 1))
#     # so output (in value terms) minus stuff being purchased by others (value terms so trade costs)
#     # per line ~ 70 below, if we sum down a row this is the world demand of a countries commodity. 

#     asset_market = A_demand

#     ##################################################################
#     # Run gravity regression on model "data"

#     trademodel = log.(normalize_by_home_trade(πshare, Ncntry)')

#     dfmodel = hcat(DataFrame(trade = vec(drop_diagonal(trademodel, Ncntry))), grv_params.dfcntryfix)

#     grvmodel = gravity(dfmodel, trade_cost_type =  trade_cost_type)

#     out_moment_vec = [grvmodel.S[1:end-1] .- grvdata.S[1:end-1] ; 
#         grvmodel.θm[1:end-1] .- grvdata.θm[1:end-1] ;
#         grvmodel.dist_coef .- grvdata.dist_coef;
#         grvmodel.lang_coef .- grvdata.lang_coef]

#     ##################################################################

#     return [sum(A_demand); goods_market[2:end]; out_moment_vec]

# end

##################################################################