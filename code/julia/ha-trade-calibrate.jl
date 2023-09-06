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

    πshare, hh, dist = world_equillibrium(Rsol, Wsol, foo_hh_params, calibrate_cntry_params)[5:7];
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

    return out_moment_vec, Wsol, β, πshare, hh, dist
end

# #####################################################################################################

function calibrate_micro(xxx, R, micro_moment, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")
    # multiple dispatch version that creates own initial guess

    dfWguess = DataFrame(CSV.File("wage-guess.csv"))

    initial_x = log.([dfWguess.guess[1:18]; 0.93])

    gravx = xxx[1:end-1]

    foo_hh_params = household_params(hh_params, ψslope = xxx[end])

    out_macro_vec, Wsol, β, πshare, hh, dist = calibrate(gravx, R, initial_x, grvdata, grvparams, foo_hh_params, cntry_params; 
            tol_vfi = tol_vfi, tol_dis = tol_dis, 
            hh_solution_method = hh_solution_method, 
            stdist_sol_method = stdist_sol_method , trade_cost_type = trade_cost_type)

    new_hh_params = household_params(hh_params, β = β)

    TFP, d = make_country_params(gravx, cntry_params, grvparams, trade_cost_type = trade_cost_type)

    #     ##################################################################
    cntry = 19 # look at the US

    p = make_p(Wsol[1:end], TFP, d[cntry, :],  cntry_params.tariff[cntry, :] )

    fooX = make_Xsection(R, Wsol[cntry], p, hh[cntry], dist[cntry], cntry,  new_hh_params; Nsims = 100000);

    micro_model = cal_πii_make_stats(fooX, prctile = [20.0, 80.0])

    out_micro_vec = [(micro_model.rich_πii - micro_model.poor_πii) - micro_moment]

    return [out_macro_vec; out_micro_vec]

end

# #####################################################################################################

# function calibrate_micro(xxx, micro_moment, R, initial_x, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
#     hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")

#     out_moment_vec, Wsol, β, πshare, hh, dist = calibrate(xxx, R, initial_x, grvdata, grvparams, hh_params, cntry_params; tol_vfi = 1e-6, tol_dis = 1e-10, 
#             hh_solution_method = hh_solution_method, stdist_sol_method = stdist_sol_method , trade_cost_type = trade_cost_type)

#     foo_hh_params = household_params(hh_params, β = β)

#     TFP, d = make_country_params(xxx, cntry_params, grvparams, trade_cost_type = trade_cost_type)

#     ##################################################################
#     cntry = 19 # look at the US

#     p = make_p(Wsol[1:end], TFP, d[cntry, :],  cntry_params.tariff[cntry, :] )

#     fooX = make_Xsection(R[cntry], Wsol[cntry], p, hh[cntry], dist[cntry], cntry, foo_hh_params; Nsims = 100000);

#     micro_model = cal_πii_make_stats(fooX, prctile = [20.0, 80.0])

#     out_micro_vec = [(micro_model.rich_πii - micro_model.poor_πii) - micro_moment]

#     return [out_moment_vec; out_micro_vec]

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

function calibrate_world_equillibrium(xxx, R, grvdata, grv_params, hh_params, cntry_params; tol_vfi = 1e-10, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")

# This is the function to try and do everything in one step, 
# still not working...

    @unpack Ncntry = cntry_params

    ## A bunch of organization here ####################

    prices = xxx[1:Ncntry]

    W = [prices[1:(Ncntry - one(Ncntry))]; 19.0 - sum(prices[1:(Ncntry - one(Ncntry))])]
    W = W./ ( sum(W / Ncntry) )

    R = ones(Ncntry)*R

    foo_hh_params = household_params(hh_params, β = xxx[Ncntry])

    TFP_grav_params = xxx[Ncntry+1:end]

    @assert length(TFP_grav_params) == 2*(Ncntry - 1) + 4 + 6  

    TFP, d = make_country_params(TFP_grav_params, cntry_params, grv_params, trade_cost_type = trade_cost_type)

    ##################################################################

    calibrate_cntry_params = country_params(TFP = TFP, d = d, Ncntry = Ncntry, L = cntry_params.L)

    Y, tradeflows, A_demand, Gbudget, πshare = world_equillibrium(R, W, foo_hh_params, calibrate_cntry_params, tol_vfi = tol_vfi, tol_dis = tol_dis, 
        hh_solution_method = hh_solution_method, stdist_sol_method=stdist_sol_method)[1:5]

    goods_market = ( Y .- vec(sum(tradeflows, dims = 1)) ) 
    # so output (in value terms) minus stuff being purchased by others (value terms so trade costs)
    # per line ~ 70 below, if we sum down a row this is the world demand of a countries commodity. 

    ##################################################################
    # Run gravity regression on model "data"

    trademodel = log.(normalize_by_home_trade(πshare, Ncntry)')

    dfmodel = hcat(DataFrame(trade = vec(drop_diagonal(trademodel, Ncntry))), grv_params.dfcntryfix)

    grvmodel = gravity(dfmodel, trade_cost_type =  trade_cost_type)

    out_moment_vec = [grvmodel.S[1:end-1] .- grvdata.S[1:end-1] ; 
        grvmodel.θm[1:end-1] .- grvdata.θm[1:end-1] ;
        grvmodel.dist_coef .- grvdata.dist_coef;
        grvmodel.lang_coef .- grvdata.lang_coef] 

    ##################################################################

    return [sum(A_demand); goods_market[2:end]; out_moment_vec]

end

##################################################################


function calibrate_world_equillibrium(xxx, R, micro_moment, grvdata, grv_params, 
    hh_params, cntry_params; tol_vfi = 1e-10, tol_dis = 1e-10, 
    hh_solution_method = "itteration", stdist_sol_method = "itteration", trade_cost_type = "ek")
    # using mulitple dispacth here
    # tries to do in all one step + micro moments (the shares)

    @unpack Ncntry = cntry_params

    ## A bunch of organization here ####################

    prices = xxx[1:Ncntry]

    W = [prices[1:(Ncntry - one(Ncntry))]; 19.0 - sum(prices[1:(Ncntry - one(Ncntry))])]
    W = W./ ( sum(W / Ncntry) )

    R = ones(Ncntry)*R

    foo_hh_params = household_params(hh_params, β = xxx[Ncntry], ψslope = xxx[Ncntry+1])

    TFP_grav_params = xxx[Ncntry+2:end]

    @assert length(TFP_grav_params) == 2*(Ncntry - 1) + 4 + 6  

    TFP, d = make_country_params(TFP_grav_params, cntry_params, grv_params, trade_cost_type = trade_cost_type)

    ##################################################################

    calibrate_cntry_params = country_params(TFP = TFP, d = d, Ncntry = Ncntry, L = cntry_params.L)

    Y, tradeflows, A_demand, Gbudget, πshare, hh, dist = world_equillibrium(R, W, foo_hh_params, calibrate_cntry_params, tol_vfi = tol_vfi, tol_dis = tol_dis, 
        hh_solution_method = hh_solution_method, stdist_sol_method=stdist_sol_method)

    goods_market = ( Y .- vec(sum(tradeflows, dims = 1)) ) 
    # so output (in value terms) minus stuff being purchased by others (value terms so trade costs)
    # per line ~ 70 below, if we sum down a row this is the world demand of a countries commodity. 

    ##################################################################
    # Run gravity regression on model "data"

    trademodel = log.(normalize_by_home_trade(πshare, Ncntry)')

    dfmodel = hcat(DataFrame(trade = vec(drop_diagonal(trademodel, Ncntry))), grv_params.dfcntryfix)

    grvmodel = gravity(dfmodel, trade_cost_type =  trade_cost_type)

    out_moment_vec = [grvmodel.S[1:end-1] .- grvdata.S[1:end-1] ; 
        grvmodel.θm[1:end-1] .- grvdata.θm[1:end-1] ;
        grvmodel.dist_coef .- grvdata.dist_coef;
        grvmodel.lang_coef .- grvdata.lang_coef] 

    ##################################################################
    cntry = 19 # look at the US

    p = make_p(W[1:end], TFP, d[cntry, :],  calibrate_cntry_params.tariff[cntry, :] )

    fooX = make_Xsection(R[cntry], W[cntry], p, hh[cntry], dist[cntry], cntry, foo_hh_params; Nsims = 100000);

    micro_model = cal_πii_make_stats(fooX, prctile = [20.0, 80.0])

    out_micro_vec = [(micro_model.rich_πii - micro_model.poor_πii) - micro_moment]

    ##################################################################

    return [sum(A_demand); goods_market[2:end]; out_moment_vec; out_micro_vec]

end