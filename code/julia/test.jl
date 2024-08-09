include("ha-trade.jl")

R = [1.014; 1.014]
T = 2
Ncntry = 2
Na = 50
Nshocks = 10

TFP = [1.0; 1.0]

τ = [0.0; 0.0]

L = [1.0; 1.0]

tariff = zeros(Ncntry, Ncntry)

d_ij = 1.745
d = [1.0 d_ij; d_ij 1.0]

Rpath = repeat(R, outer = (1,T)) 
Rend = copy(R)

R_test = reshape(Rpath, Ncntry, T)
hcat(R_test, Rend)

d_path = d.* ones(Ncntry, Ncntry, T)

cntry_prm = country_params(Ncntry = Ncntry, L = L, d = d, TFP = TFP)

Q = Array{Array{Float64,3}}(undef, T)
    #Q = Array{Float64}(undef, Na*Nshocks, Na*Nshocks)

for fwdate = 1:T

    #### Country dimension needs to be changed

    Q[fwdate] = Array{Float64}(undef, Na*Nshocks, Na*Nshocks, Ncntry)

end

size(Q)

Q[1][:,:,1]

hcat(zeros(Ncntry, T), TFP)

vcat([1.0; 2.0], [3.0; 4.0])[:]
reshape(vcat([1.0; 2.0], [3.0; 4.0]), Ncntry, T)

W_path = repeat( [1.0; 2.0; 3.0], outer = (1,T))
W_path = W_path[:]