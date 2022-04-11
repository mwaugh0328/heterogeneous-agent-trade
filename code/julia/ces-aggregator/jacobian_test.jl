function jac_test(x)
    
    outvec = Array{eltype(x)}(undef,2)
    
    outvec[1] = x[1]^2.0 * x[2]

    outvec[2] = 5 * x[1] + sin(x[2])

    return outvec

end


function exact_jac(x)

    outvec = Array{eltype(x)}(undef,2,2)

    outvec[1,1] = 2.0 * x[1] * x[2]

    outvec[1,2] = x[1] ^ 2.0

    outvec[2,1] = 5.0

    outvec[2,2] = cos(x[2])

    return outvec

end
