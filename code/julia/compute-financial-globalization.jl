include("ha-trade-environment.jl")
include("ha-trade-solution.jl")
include("ha-trade-helper-functions.jl")
include("static-trade-environment.jl")
using MINPACK
using Plots
using CSV
using DataFrames

####################################################################################
####################################################################################
println(" ")
println(" ")
println("########### computing initial eq ################")
println(" ")

Ncntry = 19

dftrade = DataFrame(CSV.File("ek-trade.csv"))

d = reshape(dftrade.d, Ncntry,Ncntry)

df = DataFrame(CSV.File("solution.csv"))

initial_x = [df.wage[2:end]; 1.0235]

TFP = df.TFP
L = df.L

Ncntry = size(d)[1]

mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, 
γ = 1.5, ϕ = 2.0, amax = 8.0, σϵ = 0.25, d = d, TFP = TFP, L = L)

f(x) = world_equillibrium_FG(x, mdl_prm, hh_solution_method = "itteration", stdist_sol_method = "itteration");

function f!(fvec, x)

    fvec .= f(x)

end

###################################################################

n = length(initial_x)
diag_adjust = n - 1

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

print(sol)

Wsol = [1.0; sol.x[1:(Ncntry - 1)]]

Rsol = ones(Ncntry)*sol.x[end]

Y, tradeflows, A_demand, tradeshare, hh, dist = world_equillibrium(Rsol,
    Wsol, mdl_prm, hh_solution_method = "itteration");

plot(log.(vec(tradeshare)), log.(dftrade.tradesharedata), seriestype = :scatter)

dftrade_model_data = DataFrame(
    trademodel = log.(vec(tradeshare)),
    tradedata = log.(dftrade.tradesharedata)
     );

CSV.write("../../notebooks/trade_model_data_fg.csv", dftrade_model_data)

hh_df = make_hh_dataframe(dist, hh, 19, Rsol, Wsol, mdl_prm)

CSV.write("../../notebooks/household_data_pre_fg.csv", hh_df)

# ####################################################################################
# println(" ")
# println(" ")
# println("########### computing counter factual eq ################")
# println(" ")


country = 4

Δ_d = 0.10
d_prime = deepcopy(d)
d_prime[19,country] =  (d[19,country]).*(1.0 + Δ_d)

country_name = "-"*string(country)


Δ_mdl_prm = world_model_params(Ncntry = Ncntry, Na = 100, 
γ = 1.5, ϕ = 2.0, amax = 8.0, σϵ = 0.25, d = d_prime, TFP = TFP, L = L)

# ###################################################################################
# # Fix prices, change d, see what happens...

Δp_Y, Δp_tradeflows, Δp_A_demand, Δp_tradeshare, Δp_hh, Δp_dist = world_equillibrium(Rsol,
Wsol, Δ_mdl_prm, hh_solution_method = "itteration");

∂W, ∂logW = welfare_by_state(hh, Δp_hh, 19, Δ_mdl_prm.σϵ)

dfwelfare = make_welfare_dataframe(∂W, ∂logW, Δ_mdl_prm)

root = "../../notebooks/welfare-US"
CSV.write(root*country_name*"-fix-p-fg.csv", dfwelfare)

hh_df = make_hh_dataframe(Δp_dist, Δp_hh, 19, Rsol, Wsol, Δ_mdl_prm)

root = "../../notebooks/household-data-US"
CSV.write(root*country_name*"-fix-p-fg.csv", hh_df)


# ###################################################################################
f(x) = world_equillibrium_FG(exp.(x), Δ_mdl_prm, hh_solution_method = "itteration", stdist_sol_method = "itteration");

function f!(fvec, x)

    fvec .= f(x)

end

initial_x =[ones(Ncntry-1); 1.02]

###################################################################

n = length(initial_x)
diag_adjust = n - 1

Δ_sol = fsolve(f!, log.(initial_x), show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

print(Δ_sol)

Δ_Wsol = exp.([0.0; Δ_sol.x[1:Ncntry-1]])
Δ_Rsol = ones(Ncntry)*exp.(Δ_sol.x[end])

Δ_Y, Δ_tradeflows, Δ_A_demand, Δ_tradeshare, Δ_hh, Δ_dist = world_equillibrium(Δ_Rsol,
Δ_Wsol, Δ_mdl_prm, hh_solution_method = "itteration")

∂W, ∂logW = welfare_by_state(hh, Δ_hh, 19, Δ_mdl_prm.σϵ)

dfwelfare = make_welfare_dataframe(∂W, ∂logW, Δ_mdl_prm)

root = "../../notebooks/welfare-US"
CSV.write(root*country_name*"-fg.csv", dfwelfare)

hh_df = make_hh_dataframe(Δ_dist, Δ_hh, 19, Δ_Rsol, Δ_Wsol, Δ_mdl_prm)

root = "../../notebooks/household-data-US"
CSV.write(root*country_name*"-fg.csv", hh_df)

global_trade_elasticity =  (log.(Δ_tradeshare ./ diag(Δ_tradeshare)) .- 
    log.(tradeshare ./ diag(tradeshare))) ./ (log.(d_prime) .- log.(d))

