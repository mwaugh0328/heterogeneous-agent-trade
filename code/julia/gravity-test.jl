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
# Recover the trade costs and technology parameters

θ = 4.0

d = zeros(19,19)
T = zeros(19)
W = ones(19)

@time gravity!(dftrade, d, T, W, θ)

################################################################
# Feed into the model and then compare with data

@time πshares, Φ = eaton_kortum(d, T, W, θ)

@time πdata =  make_trade_shares(dftrade)

plot(log.(vec(πshares)), log.(vec(πdata)), seriestype = :scatter)
