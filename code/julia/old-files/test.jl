importer_index = Array{Float64}(undef, Ncntry, Ncntry)
exporter_index = Array{Float64}(undef, Ncntry, Ncntry)

for importer = 1:Ncntry

    for exporter = 1:Ncntry

        importer_index[importer,exporter] = importer

        exporter_index[importer, exporter] = exporter

    end
end