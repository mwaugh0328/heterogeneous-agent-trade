struct θs{T}
    θπ::Array{T} 
    θc::Array{T} 
    θπii::Array{T} 
    θcii::Array{T}
end

##########################################################################
##########################################################################

function aggregate_θ(θs, ω, homecntry, model_params)

    @unpack θπ, θc, θπii, θcii = θs

    agθ = Array{Float64}(undef, model_params.Ncntry)

    for cntry = 1:model_params.Ncntry

        homeθ = sum(( θπii[:, :, cntry] .+ θcii[:, :, cntry] ) .* ω[:, :, homecntry])

        agθ[cntry] = 1 + sum(( θπ[:, :, cntry] .+ θc[:, :, cntry] ) .* ω[:, :, cntry]) - homeθ

    end

    return agθ

end

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

function make_θ(household, homecontry, R, W, p, model_params; points = 3, order = 1)

    @unpack cons_policy, Tv = household

    θπ = similar(cons_policy)
    θc = similar(cons_policy)
    θπii = similar(cons_policy)
    θcii = similar(cons_policy)

    @inbounds for idxj = 1:model_params.Ncntry

        # first construct the intensive margin elascitity

        h(x) = make_θπ(x, homecontry, idxj, p, cons_policy, Tv, R, W, model_params)
        
        foobar = central_fdm(points, order, adapt = 0)(h, log.(p[idxj]))

        θπ[:,:,idxj] .= foobar[1:model_params.Na,:] 

        θπii[:,:,idxj] .= foobar[model_params.Na+1:end,:] 

    end

    @inbounds for idxj = 1:model_params.Ncntry

        # first construct the intensive margin elascitity

        h(x) = make_θc(x, homecontry, idxj, p, cons_policy, Tv, R, W, model_params)
        
        foobar = central_fdm(points, order, adapt = 0)(h, log.(p[idxj]))

        θc[:,:,idxj] .= foobar[1:model_params.Na,:] 

        θcii[:,:,idxj] .= foobar[model_params.Na+1:end,:] 

    end

    return θs(θπ, θc, θπii, θcii)

end

##########################################################################
##########################################################################


function make_θπ(logp, homecntry, idxj, pvec, gc, v, R, W, model_params)
    # function to compute extensive margin elasticities 

    Tv = similar(v)
    
    foo = copy(pvec)
    
    foo[idxj] = exp.(logp) #p is assumed to be in log, convert to levels
    
    Tv .= coleman_operator(gc, v, R, W, foo, model_params)[2]
    # find how value function changes. Note that the way this is computed
    # shares are at old ones, so this is as if only change is 
    # (i) current utitlity
    # (ii) and ∂V/∂a (which should be zero via envelope theorem)
    # no change from future shift in shares

    #∂shares then computes **only** numerator of share, so the assumption
    #here is denominator is fixed, this is mc assumption

    πprob = make_πprob(Tv, model_params.σϵ)
    
    return vcat(log.( πprob[:,:,idxj] ), log.( πprob[:,:,homecntry] ))
    
end

##########################################################################
##########################################################################

function make_θc(logp, homecntry, idxj, pvec, gc, v, R, W, model_params)
    # function to compute intensive margin elasticities numerically

    Kgc = similar(gc)
    
    foo = copy(pvec)
    
    foo[idxj] = exp.(logp) #p is assumed to be in log, convert to levels
    
    Kgc .= coleman_operator(gc, v, R, W, foo, model_params)[1]
    
    return vcat(log.( Kgc[:,:,idxj] ), log.( Kgc[:,:,homecntry] ))
    
end