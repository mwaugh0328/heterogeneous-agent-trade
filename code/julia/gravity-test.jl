include("gravity-tools.jl")
include("trade-environment.jl")
using CSV
using DataFrames
using Plots

################################################################
# builds the EK dataset

dftrade = DataFrame(CSV.File("../../ek-data/ek-data.csv"))

dflang = DataFrame(CSV.File("../../ek-data/ek-language.csv"))

filter!(row -> ~(row.trade ≈ 1.0), dftrade);

filter!(row -> ~(row.trade ≈ 0.0), dftrade);

dftrade = hcat(dftrade, dflang);

################################################################
# Run the Gravity regression

gravity(dftrade, display = true);

################################################################
# Compute eq. wages from trade imbalance 

πdata =  make_trade_shares(dftrade)

f(x) = trade_balance(x, dflabor.L, πdata)

function f!(fvec, x)

    fvec .= f(x)

end

initial_x = ones(19)

n = length(initial_x)
diag_adjust = n - 1

sol = fsolve(f!, initial_x, show_trace = true, method = :hybr;
      ml=diag_adjust, mu=diag_adjust,
      diag=ones(n),
      mode= 1,
      tol=1e-5,
       )

print(sol)

W = sol.x ./ sol.x[19]


################################################################
# Recover the trade costs and technology parameters

θ = 4.0

d = zeros(19,19)
T = zeros(19)


@time gravity!(dftrade, d, T, W, θ)

################################################################
# Feed into the model and then compare with data

@time πshares, Φ = eaton_kortum(W, d, T, θ)

πdata =  make_trade_shares(dftrade)

plot(log.(vec(πshares)), log.(vec(πdata)), seriestype = :scatter)
