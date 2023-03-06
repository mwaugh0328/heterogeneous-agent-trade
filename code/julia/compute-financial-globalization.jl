include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("ha-trade-elasticity.jl")
include("ha-trade-efficient.jl")
include("static-trade-environment.jl")
include("gravity-tools.jl")
using MINPACK
using Plots
using CSV
using DataFrames

####################################################################################
####################################################################################
# This sets up the EK trade data and gravity stuff

dftrade = DataFrame(CSV.File("../../ek-data/ek-data.csv"))

dftrade.trade = parse.(Float64, dftrade.trade)
    # forsome reason, now it reads in as a "String7"
    
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

df = DataFrame(CSV.File("solution-fg.csv"))

L = df.L

Ncntry = size(L)[1]

γ = 1.5
σϵ = 0.25

hh_prm = household_params(Ncntry = Ncntry, Na = 100, β = 0.92,
γ = γ, ϕ = 0.5, amax = 8.0, σϵ = σϵ)

cntry_prm = country_params(Ncntry = Ncntry, L = L)

dfparams = DataFrame(CSV.File("current-guess-ek-scale.csv"))
#dfparams = DataFrame(CSV.File("current-guess-log-ek.csv"))

xxx = dfparams.guess

out_moment_vec, Wsol, Rsol, πshare  = calibrate(xxx, grvdata, grv_params, hh_prm, cntry_prm, trade_cost_type = trade_cost_type)

TFP, d = make_country_params(xxx, cntry_prm, grv_params, trade_cost_type = trade_cost_type)

cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP)

####################################################################################
####################################################################################


Y, tradeflows, A_demand, tradeshare, hh, dist = world_equillibrium(Rsol, Wsol, hh_prm, cntry_prm);

# This is a Plot test to make sure this is doing what I think it is

trademodel = log.(vec(normalize_by_home_trade(tradeshare, grv_params.Ncntry)'))

dfplot = DataFrame(trademodel = trademodel)

filter!(row -> ~(row.trademodel ≈ 1.0), dfplot);

filter!(row -> ~(row.trademodel ≈ 0.0), dfplot);

dfplot = hcat(dftrade, dfplot);

plot(dfplot.trademodel, dfplot.trade, seriestype = :scatter, alpha = 0.75,
    xlabel = "model",
    ylabel = "data",
    legend = false)

rootfile = "../../notebooks/output/"
day= "0306.csv"

home_country = 19

hh_df = make_hh_dataframe(dist, hh, home_country, Rsol, Wsol, hh_prm)

CSV.write(rootfile*"household_data_pre_fg"*day, hh_df)

# # ####################################################################################
println(" ")
println(" ")
println("########### computing counter factual eq ################")
println(" ")

country = 8
country_name = "-JPN"

Δ_d = 0.10

d_prime = deepcopy(d)
d_prime[home_country, country] =  (d[home_country, country]).*(1.0 - Δ_d)


Δ_cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d_prime, TFP = TFP)

###################################################################################
# Fix prices, change d, see what happens...

Δp_Y, Δp_tradeflows, Δp_A_demand, Δp_tradeshare, Δp_hh, Δp_dist = world_equillibrium(Rsol, Wsol, 
            hh_prm, Δ_cntry_prm, tol_vfi = 1e-10);

∂W, ∂logW = welfare_by_state(hh, Δp_hh, home_country, hh_prm.σϵ*(TFP[home_country]^(1.0 - γ)))

dfwelfare = make_welfare_dataframe(∂W, ∂logW, hh_prm)

# root = rootfile*"welfare-US"

# CSV.write(root*country_name*"-fix-p-fg"*day, dfwelfare)

# hh_df = make_hh_dataframe(Δp_dist, Δp_hh, home_country, Rsol, Wsol, hh_prm)

# root = rootfile*"household-data-US"

# CSV.write(root*country_name*"-fix-p-fg"*day, hh_df)


# # ###################################################################################

f(x) = world_equillibrium_FG(exp.(x), hh_prm, Δ_cntry_prm)

function f!(fvec, x)

    fvec .= f(x)

end

# ###################################################################

xguess = [Wsol[1:18]; Rsol[1]]

n = length(xguess)
diag_adjust = n - 1

Δ_sol = fsolve(f!, log.(initial_x), show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

print(Δ_sol)

Δ_Wsol = exp.([Δ_sol.x[1:Ncntry-1] ; 0.0 ])
Δ_Rsol = ones(Ncntry)*exp.(Δ_sol.x[end])

Δ_Y, Δ_tradeflows, Δ_A_demand, Δ_tradeshare, Δ_hh, Δ_dist = world_equillibrium(Δ_Rsol,
Δ_Wsol, hh_prm, Δ_cntry_prm, tol_vfi = 1e-10);

∂W, ∂logW = welfare_by_state(hh, Δ_hh, home_country, hh_prm.σϵ*(TFP[home_country]^(1.0 - γ)))

dfwelfare = make_welfare_dataframe(∂W, ∂logW, hh_prm)

# root = rootfile*"welfare-US"
# CSV.write(root*country_name*"-ge-fg-log.csv", dfwelfare)

# hh_df = make_hh_dataframe(Δ_dist, Δ_hh, 19, Δ_Rsol, Δ_Wsol, Δ_mdl_prm)

# root = rootfile*"household-data-US"
# CSV.write(root*country_name*"-ge-fg-log.csv", hh_df)

# global_trade_elasticity =  (log.(Δ_tradeshare ./ diag(Δ_tradeshare)) .- 
#     log.(tradeshare ./ diag(tradeshare))) ./ (log.(d_prime) .- log.(d))

