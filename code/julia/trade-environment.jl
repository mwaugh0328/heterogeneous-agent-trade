
function eaton_kortum(d, T, W, θ)
    # constructs pattern of trade for eaton and kortum model

    Ncntry = size(d)[1]
    
    πshares = Array{eltype(d)}(undef, size(d))
    
    Φ = similar(T)

    for importer = 1:Ncntry

        for exporter = 1:Ncntry

            πshares[importer, exporter] = T[exporter] * (W[exporter] * d[importer, exporter])^(-θ)
            # equation (8)-like 

        end

        Φ[importer] = sum( πshares[importer, :])
        # equation (7)

        πshares[importer, : ] .= πshares[importer, : ] / Φ[importer]
        # complete equation (8)

    end

    return πshares, Φ

end

