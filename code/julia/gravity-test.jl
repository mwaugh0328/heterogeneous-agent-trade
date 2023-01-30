using CSV
using DataFrames

include("gravity-tools.jl")

dftrade = DataFrame(CSV.File("../../ek-data/ek-data.csv"))

dflang = DataFrame(CSV.File("../../ek-data/ek-language.csv"))

filter!(row -> ~(row.trade ≈ 1.0), dftrade);

filter!(row -> ~(row.trade ≈ 0.0), dftrade);

dftrade = hcat(dftrade, dflang);