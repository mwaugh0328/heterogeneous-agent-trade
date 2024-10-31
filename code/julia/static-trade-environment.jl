function logit_trade(W, p, model_params)

    @unpack mc, γ, σϵ, Na, Nshocks, Ncntry = model_params

    shocks = exp.(mc.state_values)

    gc = Array{eltype(R)}(undef, Na, Nshocks, Ncntry)
    v = similar(gc)

    @views @inbounds for cntry = 1:Ncntry
        # fix the country
        for shk = 1:Nshocks

            for ast = 1:Na

                gc[ast,shk,cntry] = (W / p[cntry])*shocks[shk]

                v[ast,shk,cntry] = utility_fast(gc[ast, shk, cntry], γ)

            end

        end

    end

    πprob = make_πprob(v, σϵ)

    return gc, πprob, v

end

function basic_CES(W, p, σ, model_params)

    @unpack mc, γ, σϵ, Na, Nshocks, Ncntry = model_params

    shocks = exp.(mc.state_values)

    gc = Array{eltype(W)}(undef, Na, Nshocks, Ncntry)
    v = similar(gc)

    p_index = sum(p .^ -σ)

    @views @inbounds for cntry = 1:Ncntry
        # fix the country
        for shk = 1:Nshocks

            for ast = 1:Na

                gc[ast,shk,cntry] = (p[cntry]^(-σ-1.0) / p_index) * W * shocks[shk]

                v[ast,shk,cntry] = utility_fast(gc[ast, shk, cntry], γ)

            end

        end

    end

    πprob = make_πprob(v, σϵ)

    return gc, πprob, v

end