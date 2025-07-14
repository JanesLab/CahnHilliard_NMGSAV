using LinearAlgebra

include("laplace.jl")

function error2!(rr, c_old, c_new, mu, nxt, nyt, dt, xright, xleft, yright, yleft, boundary)
    x = 0.0

    rr .= mu .- c_old

    sor = zeros(Float64, nxt, nyt)
    laplace!(sor, rr, nxt, nyt, xright, xleft, yright, yleft, boundary)

    rr .= sor .- (c_new .- c_old) ./ dt

    for i in 1:nxt
        for j in 1:nyt
            x += rr[i, j]^2
        end
    end
    res2 = sqrt(x / (nxt * nyt))
    return res2
end
