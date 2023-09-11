include("ha-trade.jl")
include("ha-trade-welfare.jl")
using MINPACK
using Plots
using CSV
using DataFrames
using StatsBase

####################################################################################
####################################################################################
# This sets up the EK trade data and gravity stuff

cntryname = ["AUS","AUT","BEL","CAN","DNK","FIN","FRA",
    "DEU","GRC","ITA", "JPN", "NLD", "NZL","NOR","PRT",
    "ESP","SWE","GBR","USA"]

dftrade = DataFrame(CSV.File("../../ek-data/ek-data.csv"))

dftrade.trade = parse.(Float64, dftrade.trade)
    # for some reason, now it reads in as a "String7"
    
dflang = DataFrame(CSV.File("../../ek-data/ek-language.csv"))
    
dflabor = DataFrame(CSV.File("../../ek-data/ek-labor.csv"))
    
filter!(row -> ~(row.trade ≈ 1.0), dftrade);
    
filter!(row -> ~(row.trade ≈ 0.0), dftrade);
    
dftrade = hcat(dftrade, dflang);
    
    #dfcntryfix = select(dftrade,Not("trade"))
dfcntryfix = DataFrame(CSV.File("../../ek-data/ek-cntryfix.csv"))
    # these are the fixed characteristics of each country...


trade_cost_type = "ek"

grvdata = gravity(dftrade, display = true, trade_cost_type = trade_cost_type );

trc = trade_costs(grvdata.dist_coef, grvdata.lang_coef, grvdata.θm)

grv_params = gravity_params(L = dflabor.L, dfcntryfix = dfcntryfix, Ncntry = 19)

####################################################################################
####################################################################################
# Compute the EQ at the gravity parameters

dfparams = DataFrame(CSV.File("./calibration-files/current-guess-145-30-725.csv"))

xxx = dfparams.guess

L = dflabor.L

Ncntry = size(L)[1]

γ = 1.450
σϵ = 0.33
ψslope = 0.725

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ, ψslope = ψslope)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

R = 1.01

out_moment_vec, Wsol, β, πshare  = calibrate(xxx, R, grvdata, grv_params, hh_prm, 
                            cntry_prm, trade_cost_type = trade_cost_type)

TFP, d = make_country_params(xxx, cntry_prm, grv_params, trade_cost_type = trade_cost_type)

cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP)

hh_prm = household_params(hh_prm, β = β)

Rsol = ones(Ncntry)*R

# ####################################################################################
# ####################################################################################

Y, tradeflows, A_demand, Gbudget, tradeshare, hh, dist = world_equillibrium(Rsol, Wsol, hh_prm, cntry_prm, tol_vfi = 1e-10);

trademodel = log.(vec(normalize_by_home_trade(tradeshare, grv_params.Ncntry)'))

dfplot = DataFrame(trademodel = trademodel)

filter!(row -> ~(row.trademodel ≈ 1.0), dfplot);

filter!(row -> ~(row.trademodel ≈ 0.0), dfplot);

dfplot = hcat(dftrade, dfplot);

plot(dfplot.trademodel, dfplot.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)

home_country = 19

hh_df = make_hh_dataframe(dist, hh, home_country, Rsol, Wsol, hh_prm)

####################################################################################
####################################################################################
println(" ")
println(" ")
println("########### computing counter factual eq ################")
println(" ")

for lib_cntry = 9:18

    country_name = cntryname[lib_cntry]


    Δ_d = 0.10

    d_prime = deepcopy(d)
    d_prime[home_country, lib_cntry] =  (d[home_country, lib_cntry]).*(1.0 - Δ_d)
    d_prime[home_country, home_country] = 1.0


    Δ_cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d_prime, TFP = TFP)

    f(x) = world_equillibrium_FG(exp.(x), hh_prm, Δ_cntry_prm)

    function f!(fvec, x)

        fvec .= f(x)

    end


    xguess = [Wsol[1:18]; Rsol[1]]

    n = length(xguess)
    diag_adjust = n - 1

    Δ_sol = fsolve(f!, log.(xguess), show_trace = true, method = :hybr;
        ml=diag_adjust, mu=diag_adjust,
        diag=ones(n),
        mode= 1,
        tol=1e-10,
       )

    print(Δ_sol)

    foobar = Δ_sol.x    
    Δ_Wsol, Δ_Rsol = exp.([foobar[1:Ncntry-1] ; 0.0 ]) , ones(Ncntry)*exp.(foobar[end])

    Δ_Wsol = Δ_Wsol ./ ( sum(Δ_Wsol / Ncntry) )

    writedlm("./output/Δ_Wsol-"*cntryname[lib_cntry]*".txt", Δ_Wsol)
    writedlm("./output/Δ_Rsol-"*cntryname[lib_cntry]*".txt", Δ_Rsol)

    Δ_Y, Δ_tradeflows, Δ_A_demand, Δ_Gbudget, Δ_tradeshare, Δ_hh, Δ_dist = world_equillibrium(Δ_Rsol, Δ_Wsol, 
            hh_prm, Δ_cntry_prm, tol_vfi = 1e-10);


    ACR = 100*(1.0 / 4.22)*log(tradeshare[home_country,home_country] / Δ_tradeshare[home_country,home_country] )

    println(" ")
    println(" ")
    println("ACR-gains")
    println(ACR)

####################################################################################
####################################################################################
# now construct welfare and micro-moments

    ψ = make_ψ(home_country, ψslope.*TFP[home_country].^(1.0 - γ), hh_prm)
                        
    agrid = make_agrid(hh_prm, TFP[home_country])
                        
    foo_hh_prm = household_params(hh_prm, agrid = agrid, 
                TFP = TFP[home_country], L = L[home_country], σϵ = σϵ*(TFP[home_country]^(1.0 - γ)), ψ = ψ)
            

    # create **old** prices from perspective of home country
    p = make_p(Wsol[1:end], TFP, d[home_country, :], cntry_prm.tariff[home_country, :] )
            
    R = Rsol[home_country]
            
    W = Wsol[home_country]

    # construct welfare, porportional increase in total income 
    # needed at the **old** prices to match **new** value function            
    λτeqv =  eq_variation_porportional(R, W, p, Δ_hh[home_country], dist[home_country].state_index, foo_hh_prm)[1]

    writedlm("./output/welfare-"*cntryname[lib_cntry]*".txt", λτeqv)

    τsol = zeros(Δ_cntry_prm.Ncntry)

    # # compute elasticities and mpcs (can do at old prices)
    θ = make_θ(home_country, R, W, p, 
        τsol[home_country], foo_hh_prm; points = 3, order = 1)

    mpc = make_mpc(hh[home_country], R, W, p, 0.016/2, foo_hh_prm)

# # do this at the old stuff...so everything is consistent
    fooX = make_Xsection(R, W, p, hh[home_country], dist[home_country],
           θ, mpc, λτeqv, home_country, foo_hh_prm; Nsims = 100000)

    df = DataFrame(income = fooX.income, assets = fooX.a,
            homeshare = fooX.homeshare, expenditure = fooX.pc,
            mpc = fooX.mpc_avg, θ = fooX.θavg, ∂W = fooX.welfare);

    rootfile = "../../notebooks/output/"
 
    root = rootfile*"welfare-cross-section-"*cntryname[lib_cntry]*".csv"

    CSV.write(root, df);

end








 