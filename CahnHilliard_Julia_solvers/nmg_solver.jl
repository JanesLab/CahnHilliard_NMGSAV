using LinearAlgebra
using StaticArrays

include("relax.jl")
include("laplace.jl")
include("ch_error2.jl")

function source(c_old, nx, ny, dt, xright, xleft, yright, yleft, boundary)
    src_mu = zeros(Float64, nx, ny)
    src_c = zeros(Float64, nx, ny)
    ct = zeros(Float64, nx, ny)
    laplace!(ct, c_old, nx, ny, xright, xleft, yright, yleft, boundary)

    src_c .= c_old ./ dt .- ct

    return src_c, src_mu
end


function restrict_ch(uf, vf, nxc, nyc)
    uc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    vc = zeros(Float64, round(Int64, nxc), round(Int64, nyc))
    for i in 1:nxc
        for j in 1:nyc
            uc[i, j] = 0.25 * (uf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i - 1), round(Int, 2 * j)] + uf[round(Int, 2 * i), round(Int, 2 * j - 1)] + uf[round(Int, 2 * i), round(Int, 2 * j)])
            vc[i, j] = 0.25 * (vf[round(Int, 2 * i - 1), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i - 1), round(Int, 2 * j)] + vf[round(Int, 2 * i), round(Int, 2 * j - 1)] + vf[round(Int, 2 * i), round(Int, 2 * j)])
        end
    end
    return uc, vc
end

function nonL!(lap_c, lap_mu, c_new, mu_new, nxt, nyt, dt, epsilon2, xright, xleft, yright, yleft, boundary)
    ru = zeros(Float64, nxt, nyt)
    rw = zeros(Float64, nxt, nyt)
    laplace!(lap_c, c_new, nxt, nyt, xright, xleft, yright, yleft, boundary)
    laplace!(lap_mu, mu_new, nxt, nyt, xright, xleft, yright, yleft, boundary)

    for i in 1:nxt
        for j in 1:nyt
            ru[i, j] = c_new[i, j] / dt - lap_mu[i, j]
            rw[i, j] = mu_new[i, j] / dt - c_new[i, j]^3 + epsilon2 * lap_c[i, j]
        end
    end
    return ru, rw
end

# df(c) function: c.^3
# d2f(c) function: 3*c.^2
function defect!(lap_c, lap_mu, uf_new, wf_new, suf, swf, nxf, nyf, uc_new, wc_new, nxc, nyc, dt, epsilon2, xright, xleft, yright, yleft, boundary)
    ruc, rwc = nonL!(lap_c, lap_mu, uc_new, wc_new, nxc, nyc, dt, epsilon2, xright, xleft, yright, yleft, boundary)
    ruf, rwf = nonL!(lap_c, lap_mu, uf_new, wf_new, nxf, nyf, dt, epsilon2, xright, xleft, yright, yleft, boundary)
    ruf = suf - ruf
    rwf = swf - rwf
    rruf, rrwf = restrict_ch(ruf, rwf, nxc, nyc)
    duc = ruc + rruf
    dwc = rwc + rrwf
    return duc, dwc
end

function prolong_ch(uc, vc, nxc, nyc)
    uf = zeros(Float64, 2 * nxc, 2 * nyc)
    vf = zeros(Float64, 2 * nxc, 2 * nyc)
    for i in 1:nxc
        for j in 1:nyc
            uf[2*i-1, 2*j-1] = uc[i, j]
            uf[2*i-1, 2*j] = uc[i, j]
            uf[2*i, 2*j-1] = uc[i, j]
            uf[2*i, 2*j] = uc[i, j]
            vf[2*i-1, 2*j-1] = vc[i, j]
            vf[2*i-1, 2*j] = vc[i, j]
            vf[2*i, 2*j-1] = vc[i, j]
            vf[2*i, 2*j] = vc[i, j]
        end
    end
    return uf, vf
end

function vcycle!(uf_new, wf_new, su, sw, nxf, nyf, ilevel, c_relax, xright, xleft, yright, yleft, dt, epsilon2, n_level, boundary)
    relax!(uf_new, wf_new, su, sw, nxf, nyf, c_relax, xright, xleft, yright, yleft, dt, epsilon2, boundary)

    if ilevel < n_level
        nxc = trunc(Int64, nxf / 2)
        nyc = trunc(Int64, nyf / 2)
        uc_new, wc_new = restrict_ch(uf_new, wf_new, nxc, nyc)
        lap_c = zeros(Float64, nxf, nyf)
        lap_mu = zeros(Float64, nxf, nyf)
        duc, dwc = defect!(lap_c, lap_mu, uf_new, wf_new, su, sw, nxf, nyf, uc_new, wc_new, nxc, nyc, dt, epsilon2, xright, xleft, yright, yleft, boundary)

        uc_def = copy(uc_new)
        wc_def = copy(wc_new)

        # uc_def, wc_def = 
        vcycle!(uc_def, wc_def, duc, dwc, nxc, nyc, ilevel + 1, c_relax, xright, xleft, yright, yleft, dt, epsilon2, n_level, boundary)


        uc_def = uc_def - uc_new
        wc_def = wc_def - wc_new

        uf_def, wf_def = prolong_ch(uc_def, wc_def, nxc, nyc)

        uf_new = uf_new + uf_def
        wf_new = wf_new + wf_def

        relax!(uf_new, wf_new, su, sw, nxf, nyf, c_relax, xright, xleft, yright, yleft, dt, epsilon2, boundary)

    end
    # return uf_new, wf_new
end


function nmg_solver!(rr, c_old, c_new, mu, nx, ny, dt, solver_iter, tol, c_relax, xright, xleft, yright, yleft, epsilon2, n_level, boundary; suffix="", printres=true)
    it_mg2 = 0
    resid2 = 1
    sc, smu = source(c_old, nx, ny, dt, xright, xleft, yright, yleft, boundary)

    while resid2 > tol && it_mg2 < solver_iter
        # c_new, mu = 
        vcycle!(c_new, mu, sc, smu, nx, ny, 1, c_relax, xright, xleft, yright, yleft, dt, epsilon2, n_level, boundary)

        resid2 = error2!(rr, c_old, c_new, mu, nx, ny, dt, xright, xleft, yright, yleft, boundary)
        if printres
            print_mat("$(outdir)/residual_$(nx)_$(tol)_$(suffix).txt", resid2)
        end
        it_mg2 += 1
    end

    return c_new
end