struct planner{T}
    rc::Array{T} # asset_policy
    c::Array{T} # asset_policy
    Y::Array{T} # choice probabilities
    πprob::Array{T} # value function
    imports::Array{T} # value function
    tradeshare::Array{T} # value function
    expenditure::Array{T} # value function
end

##########################################################################

function efficient_equillibrium(x, model_params)
    # just a wraper to put stuff into solver

    social_planner =  compute_efficient(x, model_params)

    return social_planner.rc

end

##########################################################################

function compute_efficient(cii, model_params)

    @unpack mc, Ncntry, Nshocks, TFP, d, γ, σϵ = model_params

    χ = similar(cii)
    Y = similar(cii)
    rc = similar(cii)
    sales = similar(cii)
    expenditure = similar(cii)

    c = Array{Float64}(undef,Ncntry,Ncntry)
    πprob = Array{Float64}(undef,Ncntry,Ncntry)
    imports = similar(πprob)
    tradeshare = similar(πprob)


    ### first grab the country specific multipliers 
    # from first order condition

    N = mean_z(model_params) # number of effecieny units

    for cntry = 1:Ncntry

        χ[cntry] = muc(cii[cntry], γ) #multiplier

        Y[cntry] = TFP[cntry] * L[cntry] * N #production

    end

    # then recover the bilateral levels of consumption

    for impr = 1:Ncntry

        for expr = 1:Ncntry

            c[impr, expr] = muc_inverse( χ[expr] * d[impr, expr] , γ)

        end

            πprob[impr, :]  = make_πprob_efficient( c[impr, : ] , σϵ, γ)
    
    end

    # then construct the resource constraint

    for expr = 1:Ncntry

        sales[expr] = sum( c[ :, expr].*πprob[ :, expr].*d[ :, expr].*L[expr] ) 
        # across all importing countries, compute how much they are eating, this is c*pi*d*L
        # need to include trade costs as that's how much it takes to deliver c

        rc[expr] = Y[expr] - sales[expr] 
        # then the resource constraint is just Y - sales
    end

    # now construct the bilater trade matrix

    for impr = 1:Ncntry

        for expr = 1:Ncntry

            imports[impr, expr] = (χ[expr] * d[impr, expr] 
                                        * c[impr, expr] * πprob[impr, expr] * L[impr])

            #in valueterms, so χ_{jj}*d_{ij} is the shadow price of a unit of goods
            #j in country i

        end

        expenditure[impr] = sum(imports[impr, :])

        tradeshare[impr, :] = imports[impr, :] ./ expenditure[impr]

    end


    return planner(rc, c, Y, πprob, imports, tradeshare, expenditure)
    
end

##########################################################################

function make_πprob_efficient(c, σ, γ)

    foo = ( utility.(c, γ) - muc.(c, γ).*c )

    @fastmath foo .= @. exp( foo / σ )

    return foo ./ sum( foo ) 
   
end