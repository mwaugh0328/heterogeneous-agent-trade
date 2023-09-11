#include("armington-trade-environment.jl")
#include("armington-trade-solution.jl")
#include("durables_environment.jl")
#include("durables_solution.jl")
#include("durables_utils.jl")

R = 1.01

W = [1.0 ,1.0, 1.0]

τ_revenue = [0.0 , 0.0, 0.0]

τ = [0.0  0.0 0.2 ; 0.0  0.0 0.2 ; 0.2  0.2  0.0] # tariff rate
d = [1.0  1.0 1.5; 1.0 1.0 1.5; 1.25 1.75 1.0] # trade cost
A = [1.0, 1.0, 1.0]

trp = trade_params(A = A, τ = τ, d = d)

outvec = fullmodel_eq(W, τ_revenue, R, model_params(), trp)

println(" ")
println(outvec)
println(" ")

Ncntry = 3

function f!(fvec, x)
    
    fvec[:] = fullmodel_eq(x[1:Ncntry],
            x[Ncntry+1:end], R, model_params(), trp)

end

initial_x = [0.95 , 0.95, 0.95, 0.0, 0.0, 0.0]
n = length(initial_x)
diag_adjust = n - 1

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

println(" ")
println(sol)
println(" ")

outvec, Pces, AD, N, LFP = fullmodel_eq(sol.x[1:Ncntry], 
    sol.x[Ncntry+1:end], R, model_params(), trp, output = "all");


trade_value, trade_share, τ_revenue,
    value_demand, Pindex = compute_trade_eq(sol.x[1:Ncntry],
             AD, N, sol.x[Ncntry+1:end], trp; output = "all")