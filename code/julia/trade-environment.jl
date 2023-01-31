
function eaton_kortum(W, d, T, θ)
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

################################################################
################################################################

function trade_balance(W, L, πshares; method = "solver")
    # function to compute trade (im)balance
    trade_balance = similar(L)
    Ncntry = length(L)

    for importer = 1:Ncntry

        trade_balance[importer] = W[importer]*L[importer] - sum(πshares[:, importer] .* W .* L)
        # equation (20) in Eaton and Kortum (no exogenous income, β = 1)
        # need better explanation of this

    end

    if method == "solver"

        return trade_balance

    else

        return trade_balance

    end


end


