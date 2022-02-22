using BenchmarkTools, SpecialFunctions
using Statistics
using Parameters

# trade_params = @with_kw (
#     θ = 4.0,
#     τ = [0.0  0.0; 0.0 0.0], # tariff rate
#     d = [1.0  1.5; 1.5 1.0], # trade cost
#     A = [1.0, 1.0], #TFP
#     Ncntry = length(A),
    
# )

@with_kw struct trade_params
    θ::Float64 = 4.0
    τ::Array{Float64} = [0.0  0.0; 0.0 0.0]
    d::Array{Float64} = [1.0  1.5; 1.5 1.0]
    A::Array{Float64} = [1.0, 1.0]
    Ncntry::Int64 = length(A)
end

# simple parameter settings to test stuff out on....

##########################################################################
##########################################################################

function goods_prices(w, trade_params)
    # constructs a matrix of prices for each countries 
    # commodity given tariffs, trade costs, and unit costs
    # returns the prices and the price index
    # prices are normalized so Pindex in first country = 1

    @unpack τ, d, A, Ncntry, θ = trade_params
    @assert length(w) == Ncntry
    
    p = Array{eltype(w)}(undef, Ncntry, Ncntry)
    Pindex = Array{eltype(w)}(undef, Ncntry)

    for buyr = 1:Ncntry # buyer

        for suplr = 1:Ncntry #suppleir

            p[buyr, suplr] = (1.0 + τ[buyr, suplr]) * d[buyr, suplr] * (w[suplr] / A[suplr])

        end

    end

    p, Pindex = ces(p, θ)
    # then this renormalizes everything

    return p, Pindex

end


##########################################################################
##########################################################################
struct trade_stats{T}
    trade_value::Array{T}
    trade_share::Array{T}
    τ_revenue::Array{T}
    world_demand::Array{T}
    Pindex::Array{T}
end


##########################################################################
##########################################################################

function trade_flows(p, Pindex, AD, trade_params)
    # constructs trade flows given prices, price index, demand
    # returns flows in quantities, values, shares, tariff revenue 
    # and quantity demand of each countries commodity

    @unpack Ncntry, τ, θ = trade_params

    trade_value = Array{eltype(p)}(undef, Ncntry, Ncntry)
    trade_share = Array{eltype(p)}(undef, Ncntry, Ncntry)
    τ_revenue = Array{eltype(p)}(undef, Ncntry, Ncntry)


    for buyr = 1:Ncntry # buyer

        for suplr = 1:Ncntry #suppleir

            trade_value[buyr,suplr] = p[buyr, suplr] * ( (p[buyr, suplr] / Pindex[buyr])^(-θ) ) * AD[buyr]
            # ces demand curve multipy quantity by price to get value

            #trade_share[buyr,suplr] = trade_value[buyr,suplr] / (P[buyr] * AD[buyr])

            trade_share[buyr,suplr] = ( p[buyr, suplr] / Pindex[buyr] )^(1-θ) 

            τ_revenue[buyr,suplr] = τ[buyr, suplr] * (p[buyr, suplr] / (1.0 + τ[buyr, suplr])) * trade_value[buyr,suplr]
            # the middle adjustment is to net out of price the tariff.
            # the interpertation is that this revenue buyer country recives when importing
            # from suplier country

        end
    end

    world_demand = sum(trade_value, dims = 1)[:]
    # sum across rows, down a column gives the total amount of goods demanded from 
    # each country.

    trade = trade_stats(
        trade_value, trade_share, τ_revenue, world_demand, Pindex
        )

    return trade

end

##########################################################################
##########################################################################

function ces(p, θ)
    # function to compute CES price index.

    P = sum( p.^(1.0 - θ), dims = 2).^( 1.0 / (1.0 - θ) )
    # sum across the columns is how muc each 

    # and then normalize things so that CES price index
    # in country one is one

    p ./= P[1]

    P = sum( p.^(1.0 - θ), dims = 2).^( 1.0 / (1.0 - θ) )

    return p, P
    
end
##########################################################################
##########################################################################

function make_trade_params(d, τ, A, Ncntry)

    dmat = Array{eltype(d)}(undef,Ncntry,Ncntry)
    fill!(dmat, d)

    for xxx in 1:Ncntry
        dmat[xxx,xxx] = 1.0
    end

    τmat = Array{eltype(d)}(undef,Ncntry,Ncntry)
    fill!(τmat, τ)

    for xxx in 1:Ncntry
        τmat[xxx,xxx] = 0.0
    end

    Amat = Array{eltype(d)}(undef,Ncntry)
    fill!(Amat, A)

    return dmat, τmat, Amat
    
end