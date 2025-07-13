using FFTW

function ext(x) #good
    nx, ny = size(x)
    x_ext = zeros(2 * nx, 2 * ny)

    # Original block
    x_ext[1:nx, 1:ny] = x

    # Flip horizontally
    x_ext[1:nx, ny+1:2*ny] = x[:, end:-1:1]

    # Flip vertically
    x_ext[nx+1:2*nx, 1:ny] = x[end:-1:1, :]

    # Flip both
    x_ext[nx+1:2*nx, ny+1:2*ny] = x[end:-1:1, end:-1:1]

    return x_ext
end

function extback(x_ext) #good
    # Shrinks from 2*nx x 2*ny back to nx x ny (upper-left block)
    nx_ext, ny_ext = size(x_ext)
    nx = nx_ext ÷ 2
    ny = ny_ext ÷ 2
    x_back = x_ext[1:nx, 1:ny]

    return x_back
end


function f_SAV(phi, gamma0) #good
    # f = @(x) 0.25*(x.^2-1).^2
    fphi = (phi .^ 2 .- 1 .- gamma0) .^ 2 ./ 4
    return fphi
end

function r0_fun(phi0, hx, hy, C0, gamma0) #good
    # fphi = f_SAV(phi0, gamma0)
    E1 = fft(f_SAV(phi0, gamma0))
    tmp = imag(E1[1, 1]) == 0 ? real(E1[1, 1]) : E1[1, 1]
    r0 = sqrt(hx * hy * tmp + C0)
    return r0
end


function df(phi, gamma0) #good
    return phi .^ 3 .- (1 .+ gamma0) .* phi
end

function Lap_SAV(phi, k2, boundary)
    if boundary == "periodic"
        Lphi = real(ifft(k2 .* fft(phi)))
    elseif boundary == "neumann"
        Lphi = real(ifft(k2 .* fft_filtered(phi)))
    else
        error("Boundary condition not recognized. Use 'periodic' or 'neumann'.")
    end
    return Lphi
end

function A_inv_CN(phi, dt, k2, k4, gamma0, epsilon2, boundary) #good
    if boundary == "periodic"
        denom = 1 .+ ((dt / 2) * epsilon2 .* k4 .- (dt / 2) * gamma0 .* k2)
        Ai = real.(ifft(fft(phi) ./ denom))
    elseif boundary == "neumann"
        denom = 1 .+ ((dt / 2) * epsilon2 .* k4 .- (dt / 2) * gamma0 .* k2)
        Ai = real.(ifft(fft_filtered(phi) ./ denom))
    else
        error("Boundary condition not recognized. Use 'periodic' or 'neumann'.")
    end
    return Ai
end

function fft_filtered(x)
    return (fft(x))

    # return real(fft(x))
end

function b_fun(phi, hx, hy, C0, gamma0) #good
    E1 = fft(f_SAV(phi, gamma0))
    return df(phi, gamma0) ./ sqrt.(E1[1, 1] * hx * hy .+ C0)
end

function g_fun_CN(phi0, r0, b, dt, hx, hy, epsilon2, gamma0, C0, k2, boundary) #good
    Lap_phi0 = Lap_SAV(phi0, k2, boundary)
    Lap_Lap_phi0 = Lap_SAV(Lap_phi0, k2, boundary)

    bphi0 = fft(b .* phi0)
    bphi0 = hx * hy * bphi0[1, 1]

    E1 = fft(f_SAV(phi0, gamma0))

    g = phi0 .- (dt / 2) * epsilon2 .* Lap_Lap_phi0 .+ (dt / 2) * gamma0 .* Lap_phi0 .+
        dt .* Lap_SAV(b, k2, boundary) .* (r0 .- (1 / 4) * bphi0)

    return g
end

function r_fun(phi, phi0, r0, b, hx, hy, C0, dt, gamma0)
    bphi0 = fft(b .* phi0)
    bphi0 = hx * hy * bphi0[1, 1]

    bphi = fft(b .* phi)
    bphi = hx * hy * bphi[1, 1]

    E1 = fft(f_SAV(phi0, gamma0))

    r = r0 .+ (1 / 2) * (bphi .- bphi0)

    return r
end

function calculate_discrete_energy_sav(phi, h2, epsilon2, k2, gamma0, r, C0)

    E_gamma = h2 * fft(gamma0 / 2 * phi .^ 2)
    E_gamma = E_gamma[1, 1]

    a = r^2 - C0 - (gamma0^2 + 2 * gamma0) / 4

    a = h2 * sum(f_SAV(phi, gamma0))
    sum_i = sum((phi[2:end, :] .- phi[1:end-1, :]) .^ 2)

    b = (epsilon2 / 2) * sum_i
    sum_j = sum((phi[:, 2:end] .- phi[:, 1:end-1]) .^ 2)

    c = (epsilon2 / 2) * sum_j
    E = a + b + c + E_gamma

    return E
end

# function ch_modified_energy(phi, h2, epsilon2, k2, gamma0)
#     E1 = h2 * fft(f_SAV(phi, 0))
#     E1 = E1[1, 1]

#     E = 0.5 * epsilon2 * fft(phi .* (-real(ifft(k2 .* fft(phi)))))
#     E = E1 + E[1, 1]
#     return E
# end

function ch_r_error(r, phi, h2, C0, gamma0) #need to check if this needs to be vectorized
    E1 = fft(f_SAV(phi, gamma0))
    tmp = imag(E1[1, 1]) == 0 ? real(E1[1, 1]) : E1[1, 1]

    D = r - sqrt(h2 * tmp + C0)
    D = D / sqrt(h2 * tmp + C0)
    return D
end