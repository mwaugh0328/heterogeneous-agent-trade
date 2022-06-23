using MAT
using Plots
using CSV
using DataFrames

file = matopen("../../ek-data/ek-output.mat")
tradesharedata = read(file, "tssdmat")'
# I set my model up so row is an importer, column is exporter and 
# that accross columns should sum to one this is the oppisite so flip

d = read(file, "rtausd")'

Ncntry = 19

importer_index = Array{Float64}(undef, Ncntry, Ncntry)
exporter_index = Array{Float64}(undef, Ncntry, Ncntry)

for importer = 1:Ncntry

    for exporter = 1:Ncntry

        importer_index[importer,exporter] = importer

        exporter_index[importer, exporter] = exporter

    end
    
end

df = DataFrame(importer_index = vec(importer_index), 
    exporter_index = vec(exporter_index),
    tradesharedata = vec(tradesharedata),
    d = vec(d),
     );

CSV.write("ek-trade.csv", df)


# same deal

# L = [0.054, 0.024, 0.029, 0.094, 0.017, 0.019,
#     0.181, 0.0225, 0.025, 0.159, 0.544, 0.043, 0.010, 0.015,
#     0.026, 0.10, 0.031, 0.186, 1.0]

# TFP = [0.36,0.30,0.22,0.47,0.32,0.41,0.61,0.75,0.14,
#     0.57,0.97,0.28,0.22,0.37,0.13,0.33,0.47,0.53,1.0]

# TFP .=  TFP.^(1. / 3.6)