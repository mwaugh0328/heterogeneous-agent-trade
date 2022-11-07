function compute_efficient(cii, model_params)

    @unpack mc, Ncntry, Nshocks, TFP, d, γ = model_params

    χ = similar(cii)
    Y = similar(cii)
    c = Array{Float64}(undef,Ncntry,Ncntry)
    πprob = Array{Float64}(undef,Ncntry,Ncntry)

    ### first grab the country specific multipliers 
    # from first order condition

    N = mean_z(model_params) # number of effecieny units

    for cntry = 1:Ncntry

        χ[cntry] = utility(cii[cntry], γ) #multiplier

        Y[cntry] = TFP[cntry] * L[cntry] * N #production

    end

    # then recover the bilateral levels of consumption

    for importer = 1:Ncntry

        for exporter = 1:Ncntry

            c[importer, exporter] = muc_inverse( χ[exporter] * d[importer, exporter] , γ)

        end

            πprob[importer, :]  = make_πprob_efficient( c[importer, :] , σ, γ)
    
    end

return Y, c

end

##########################################################################

function make_πprob_efficient(c, σ, γ)

    foo = ( utility.(cii, γ) - muc.(c, γ).*c )

    @fastmath foo .= @. exp( foo / σ )

    return foo ./ sum( foo, dims = 3) 
   
end