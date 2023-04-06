struct θs{T}
    θπ::Array{T} # how πij changes
    θc::Array{T}  # how cij changes
    θπii::Array{T} # how πii changes w.r. dij
    θcii::Array{T} # how cii changes w.r. dij
end

##########################################################################
##########################################################################

function aggregate_θ(θs, ω, homecntry, model_params)

    @unpack θπ, θc, θπii, θcii = θs

    agθ = Array{Float64}(undef, model_params.Ncntry)

    for cntry = 1:model_params.Ncntry

        homeθ = sum(( θπii[:, :, cntry] .+ θcii[:, :, cntry] ) .* ω[:, :, homecntry])
        # this is the θii,j part in Proposition 3

        agθ[cntry] = 1 + sum(( θπ[:, :, cntry] .+ θc[:, :, cntry] ) .* ω[:, :, cntry]) - homeθ
        # then 1 + θ_ij - \theta_ii,j

    end

    return agθ

end

##########################################################################
##########################################################################

function make_ω(household, distribution, L, p, model_params)
    #makes the weights to aggregate elasticities

    @unpack cons_policy, πprob = household
    @unpack λ = distribution

    λ = reshape(λ, model_params.Na, model_params.Nshocks)

    ω = similar(cons_policy)

    for idxj = 1:model_params.Ncntry

        ω[:,:,idxj] .= p[idxj] .* cons_policy[:,:, idxj] .* πprob[:, :, idxj] .* λ *L

    end

    for idxj = 1:model_params.Ncntry

        ω[:,:,idxj] .= ω[:,:,idxj] / sum(ω[:,:,idxj])

    end

    return ω

end

##########################################################################
##########################################################################

function make_θ(homecontry, R, W, p, model_params; points = 3, order = 1)
    # makes the micro-level elasticities
    # using multiple-dispatch here, this is the high-level function to make 
    # everything

    θπ = Array{Float64}(undef, model_params.Na, model_params.Nshocks, model_params.Ncntry)
    θc = similar(θπ)
    θπii = similar(θπ)
    θcii = similar(θπ)

    @inbounds for idxj = 1:model_params.Ncntry

        # now construct the extensive margin elascitity

        h(x) = make_θ(x, homecontry, idxj, p, R, W, model_params)
        # this is the low level function used for differentiation
        
        foobar = jacobian(central_fdm(points, order, adapt = 0), h, log.(p[idxj]) )

        foobar = reshape(foobar[1], (model_params.Na, model_params.Nshocks, 4))

        # this then upacks things to grab both ∂πij / ∂dij and ∂πii / ∂dij
        θπ[:,:,idxj] .= foobar[:,:,1] 

        θπii[:,:,idxj] .= foobar[:,:,2] 

        θc[:,:,idxj] .= foobar[:,:,3] 

        θcii[:,:,idxj] .= foobar[:,:,4] 

    end

    return θs(θπ, θc, θπii, θcii)

end

##########################################################################
##########################################################################


function make_θ(logp, homecntry, idxj, pvec, R, W, model_params)
    # low level function to compute elasticities
    # using multiple dispatch here (see above)
    
    foo = copy(pvec)
    
    foo[idxj] = exp.(logp) #logp is assumed to be in log, convert to levels

    hh = solve_household_problem(R, W, foo, model_params, tol = 1e-10)
    # resolve the whole value function...this is unlike Mongey - Waugh, were 
    # we need to just consider a one period deviation
    
    πprob = make_πprob(hh.Tv, model_params.σϵ)
    
    return log.( πprob[:,:,idxj] ), log.( πprob[:,:,homecntry] ), 
            log.( hh.cons_policy[:,:,idxj] ), log.( hh.cons_policy[:,:,homecntry] )
    # when passed through numerical diff it returns
    # ∂πij / ∂dij, ∂πii / ∂dij, ∂cij / ∂dij, ∂cii / ∂dij
    
end
