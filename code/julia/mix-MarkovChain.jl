function mMarkovChain(nP::Int64, nIID::Int64, ρ::Float64, σP::Float64, σIID::Float64; w0=1.0)
    # This is from Manuel
    # numbers from KMP
    # ar1 = 0.9695
    # sigmaP = sqrt(0.0384)/(1.2)
    # sigmaIID = sqrt(0.0522)/(1.2)
    # P, z_vals = calibration(5, 2 , ar1, sigmaP, sigmaIID)

    # persistent component
    pP, y_vals = rouwenhorst(nP, ρ, σP, 0.0)
    y_vals = exp.(y_vals)

    # iid component
    if nIID > 1
        pIID, x_vals = rouwenhorst(nIID, 0.0, σIID, 0.0)
        x_vals = exp.(x_vals)
    else
        pIID = 1.0
        x_vals = 1.0
    end

    nZ = nP * nIID
    z_vals = Array{Float64,1}(undef, nZ)
    P = zeros(nZ, nZ)

    item = 1
    for i = 1:nP
        for j = 1:nIID
            z_vals[item] = y_vals[i] * x_vals[j]
            item += 1
        end
    end

    item = 1
    for i = 1:nP
        for j = 1:nIID
            item2 = 1
            for k = 1:nP
                for l = 1:nIID
                    P[item,item2] = pP[i,k] * pIID[j,l]
                    item2 += 1
                end
            end
            item += 1
        end
    end

    for p in eachrow(P)
        p .= p ./ sum(p)
    end

    # compute mean of ergodic
    invZ = real(inv(eigvecs(P))[nZ,:])
    invZ .= invZ ./ sum(invZ)
    meanZ = invZ' * z_vals
    # scale to get desired mean
    z_vals .= (w0 / meanZ) .* z_vals

    return MarkovChain(P, log.(z_vals))
end


# from Quantecon

# @doc doc"""
# Rouwenhorst's method to approximate AR(1) processes.
# The process follows
# ```math
#     y_t = \mu + \rho y_{t-1} + \epsilon_t
# ```
# where ``\epsilon_t \sim N (0, \sigma^2)``
# ##### Arguments
# - `N::Integer` : Number of points in markov process
# - `ρ::Real` : Persistence parameter in AR(1) process
# - `σ::Real` : Standard deviation of random component of AR(1) process
# - `μ::Real(0.0)` :  Mean of AR(1) process
# ##### Returns
# - `mc::MarkovChain{Float64}` : Markov chain holding the state values and transition matrix
# """
function rouwenhorst(N::Integer, ρ::Real, σ::Real, μ::Real=0.0)
    σ_y = σ / sqrt(1 - ρ^2)
    p  = (1 + ρ) / 2
    Θ = [p 1 - p; 1 - p p]
    ψ = sqrt(N - 1) * σ_y
    m = μ / (1 - ρ)
    state_values, p = _rouwenhorst(p, p, m, ψ, N)

    return p, state_values
end

function _rouwenhorst(p::Real, q::Real, m::Real, Δ::Real, n::Integer)
    if n == 2
        return [m - Δ, m + Δ],  [p 1 - p; 1 - q q]
    else
        _, θ_nm1 = _rouwenhorst(p, q, m, Δ, n - 1)
        θN = p * [θ_nm1 zeros(n - 1, 1); zeros(1, n)] +
             (1 - p) * [zeros(n - 1, 1) θ_nm1; zeros(1, n)] +
             q * [zeros(1, n); zeros(n - 1, 1) θ_nm1] +
             (1 - q) * [zeros(1, n); θ_nm1 zeros(n - 1, 1)]
        θN[2:end - 1, :] ./= 2
        return range(m - Δ, stop=m + Δ, length=n), θN
    end
end