
include("ha-trade.jl")
using MINPACK
using Plots
import DataFrames


haparams = model_params()
             
Ncntry = trade_params().Ncntry

new_τ = [0.00 0.10; 0.00 0.00]

T = 35

τ_path = Array{Array{Float64}}(undef, 5)
fill!(τ_path, trade_params().τ)

n_path = Array{Array{Float64}}(undef, 30)
fill!(n_path, new_τ  )

τ_path = [τ_path; n_path]

@assert length(τ_path) == T

###############################################################################################
# STEP 1 Solve for initial equillibrium

inital_tradeparams = trade_params()

f(x) = ha_trade_equilibrium(x, haparams , inital_tradeparams )

function f!(fvec, x)
    
      fvec .= f(x)
  
end

W = Array{eltype(Float64)}(undef, Ncntry)
fill!(W, 1.0)

τ_revenue = Array{eltype(Float64)}(undef, Ncntry)
fill!(τ_revenue, 0.0)

R = 1.029


initial_x = [W; τ_revenue; R]
n = length(initial_x)
diag_adjust = n - 1

solint = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

println(" ")
println(solint)
println(" ")

output_int = ha_trade_equilibrium(solint.x, haparams, inital_tradeparams, display = true)[2]

Wint, τ_rev_int, Rint = unpack_solution(solint.x, Ncntry)

dist_int = collect_intial_conditions(Wint, τ_rev_int, Rint[1], haparams, inital_tradeparams )


# ###############################################################################################
# ###############################################################################################
# # STEP 2 Solve for final ss equillibrium

tparams_end = trade_params(τ = τ_path[end])

f(x) = ha_trade_equilibrium(x, haparams , tparams_end )

function f!(fvec, x)
    
      fvec .= f(x)
  
end

sol_end = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

println(" ")
println(sol_end)
println(" ")
       
Wend, τ_rev_end, Rend = unpack_solution(sol_end.x, Ncntry)

# ###############################################################################################
# ###############################################################################################
# # STEP 3 Solve for PATH

W = repeat(Wend, T)
τ_rev = repeat(τ_rev_end, T)
R = repeat(Rend , T-1)

hh_end = collect_end_conditions(Wend, τ_rev_end, Rend[1], haparams, tparams_end )

trade_path = Array{trade_params}(undef,T)

for xxx = 1:T

      trade_path[xxx] = trade_params(τ = τ_path[xxx])

end

initial_x = [W; τ_rev; R]

f(x)  = transition_path(x, Rint, trade_path, dist_int, hh_end, haparams);

function f!(fvec, x)
    
      fvec .= f(x)
  
end

sol_path = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
       )

println(" ")
println(sol_path)
println(" ")

path_stats = transition_path(sol_path.x, Rint, trade_path, 
            dist_int, hh_end, haparams, display = true)[2]

df = make_dataset(path_stats, output_int, sol_path.x, Rint, T)








