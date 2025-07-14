
using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random
using FFTW
include("sav_solver.jl")
include("ch_initialization.jl")
include("aux_functions_SAV.jl")


Random.seed!(1234) #note that when using a random seed, you must RESTART the REPL or run in a new instance for the runs to be the same. 

"""
This function uses the nonlinear multigrid method to solve the Cahn-Hilliard equation for a specified number of time steps of size dt.

    REQUIRED INPUT
        phi0 = Initial field of chemical states in the domain, created by ch_initialization.

    OPTIONAL INPUT
        t_iter = Number of time steps simulated. Default = 1e3.
        dt = Time step. Default = 2.5e-5 characteristic times.
        solver_iter = Number of solver iterations per time step. Default = 1e4.
        tol = Solver tolerance per time step. Default = 1e-5.
        dt_out = Spacing of time steps output to phi_t as a multidimensional
           array (if less than 1e9 elements) or printed file (if greater than
           1e9 elements). Default = 1, to output every time step. 
        m = Number of mesh points over which the interface exists. Default = 4. 
        epsilon2 = Squared interfacial transition distance; if specified,
           m will be overwritten. Default = nan (do not overwrite m).
        boundary = Boundary conditions for the simulation:
           'periodic' (default) - flux on one domain border equals negative flux on the opposite border.
           'neumann' - zero flux on the domain borders.
        c_relax = Number of smoothing relaxations done at the start and end
           of each V-cycle. Default = 2.
        domain = Vector of rightmost and leftmost grid points in x and y.
           Format: [xright xleft yright yleft]. Default = [1 0 1 0].
        printres = Logical to print solver residuals to a file. Default = false.
        printphi = Logical to print phi to a file regardless of whether or 
           not it can be saved as a multidimensional array. Default = false.
         pathname = Name of the path to which phi is printed. Default = 'cd'. May include a prefix for the filename.
        C0       = Regularization parameter.
        gamma0   = Stabilization parameter.
        eta      = Relaxation parameter.
        xi_flag  = Relaxation flag (0 or 1); if 0, no relaxation (xi is set to 1).
    OUTPUT
        t_out = Time corresponding to the dt time step outputs.
        phi_t = Multidimensional array of phi over t_out.
        delta_mass_t = Vector of mass change over t_out.
        E_t = Vector of relative energy over t_out.
"""
function CahnHilliard_SAV(phi0; t_iter=1e3, dt=2.5e-5, dt_out=10, m=8, epsilon2=NaN, boundary="periodic", domain=[1 0 1 0], printres=false, printphi=false, pathname="cd", C0=1, gamma0=2, eta=0.95, xi_flag=1)
    nx, ny = size(phi0)
    xright, xleft, yright, yleft = domain
    Lx = xright - xleft
    Ly = yright - yleft

    # % Decide on the solver's mesh spacing for NEUMANN vs PERIODIC
    # %  - For Neumann: we will mirror the domain, so pass 2*hx and 2*hy into sav_solver.
    # %  - For Periodic: keep as-is.
    if boundary == "neumann"
        Lx = 2 * Lx
        Ly = 2 * Ly
        nx = 2 * nx
        ny = 2 * ny
    end
    hx = Lx / nx
    hy = Ly / ny
    h2 = hx * hy #Define mesh size


    if isnan(epsilon2)
        epsilon2 = h2 * m^2 / (2 * sqrt(2) * atanh(0.9))^2 #Define ϵ^2 if not prespecified
    else
        m = sqrt((epsilon2 * (2 * sqrt(2) * atanh(0.9))^2) / h2) #Else overwrite m
    end

    kx = 1im * vcat(0:nx÷2, -nx÷2+1:-1) * (2 * pi / Lx)
    ky = 1im * vcat(0:ny÷2, -ny÷2+1:-1) * (2 * pi / Ly)

    kxx = real(kx .^ 2)
    kyy = real(ky .^ 2)
    meshgrid(x, y) = (repeat(x, 1, length(y)), repeat(y', length(x), 1))
    kxx_mat, kyy_mat = meshgrid(kxx, kyy)

    k2 = kxx_mat .+ kyy_mat
    k4 = k2 .^ 2

    # Spectral stuff for original domain for Neumann bc to calculate energy
    # if boundary == "neumann"
    #     k_x_od = 1im * vcat(0:(nx÷2)/2, -(nx ÷ 2)/2+1:-1) * (2 * pi / (Lx / 2))
    #     k_y_od = 1im * vcat(0:(ny÷2)/2, -(ny ÷ 2)/2+1:-1) * (2 * pi / (Ly / 2))
    #     k_xx_od = real(k_x_od .^ 2)
    #     k_yy_od = real(k_y_od .^ 2)
    #     (kxx_od, kyy_od) = meshgrid(k_xx_od, k_yy_od)
    #     k2_od = real(kxx_od + kyy_od)
    # end

    # Initialization

    if boundary == "neumann"
        phi_old = ext(phi0)
    elseif boundary == "periodic"
        phi_old = copy(phi0)
    end
    phi_prev = phi_old
    r_old = r0_fun(phi_old, hx, hy, C0, gamma0) # Initialize sav state

    downsampled = nx * ny * t_iter / dt_out > 1e9 #Logical index for the need to downsample
    n_timesteps = floor(Int64, t_iter / dt_out)

    if printphi
        mass_t = zeros(Float64, n_timesteps + 1, 1)
        E_t = zeros(Float64, n_timesteps + 1, 1)
        D_t = zeros(Float64, n_timesteps + 1, 1)

        t_out = 0:dt_out*dt:(n_timesteps)*dt*dt_out
        if pathname == "cd"
            pathname = pwd()
        end
        if boundary == "neumann"
            phi_old_out = extback(phi_old)
        else
            phi_old_out = phi_old
        end

        open("$(pathname)phi.csv", "w", lock=false) do f
            writedlm(f, phi_old_out, ",")
        end
        phi_t = NaN
    else #save to variables
        if downsampled
            new_dt_out = ceil(nx * ny * t_iter / 1e9) #we need to round up to ensure we have enough space
            println("Variable phi_t is too large with dt_out = %4.0f. Downsampling to every %4.0f time steps\n", dt_out, new_dt_out)
            dt_out = new_dt_out
            n_timesteps = floor(t_iter / dt_out)
        end
        mass_t = zeros(Float64, n_timesteps + 1, 1)
        E_t = zeros(Float64, n_timesteps + 1, 1)
        D_t = zeros(Float64, n_timesteps + 1, 1)
        t_out = 0:dt_out*dt:(n_timesteps)*dt*dt_out
        if boundary == "neumann"
            phi_t = zeros(Float64, nx ÷ 2, ny ÷ 2, n_timesteps + 1)
            phi_old_out = extback(phi_old)
            phi_t[:, :, 1] = phi_old_out
        else
            phi_t = zeros(Float64, nx, ny, n_timesteps + 1)
            phi_old_out = phi_old
            phi_t[:, :, 1] = phi_old_out
        end
    end

    if boundary == "neumann"
        mass_t[1] = calculate_mass(phi0, h2, nx / 2, ny / 2)
    else
        mass_t[1] = calculate_mass(phi0, h2, nx, ny)
    end
    # E_t[1] = calculate_discrete_energy(phi_old_out, h2, epsilon2)
    # elseif boundary == "periodic"
    E_t[1] = calculate_discrete_energy(phi_old_out, h2, epsilon2)
    # end
    D_t[1] = ch_r_error(r_old, phi_old, h2, C0, gamma0)

    if printres
        println("Saving squared residuals per iteration to file in the output directory\n")
    end

    for it in 1:t_iter
        phi_new, r_new = sav_solver(phi_old, phi_prev, r_old, hx, hy, k2, k4, dt, epsilon2, boundary, C0, gamma0, eta, xi_flag, it)
        phi_new = real.(phi_new)
        r_new = real.(r_new)
        if boundary == "neumann"
            phi_new_out = extback(phi_new)
        else
            phi_new_out = phi_new
        end
        if it / t_iter * 100 % 5 == 0
            println(@sprintf("%3.0f percent complete", it / t_iter * 100))
        end
        if it % dt_out == 0
            t_index = floor(Int64, it / dt_out) + 1
            if printphi
                open("$(pathname)phi.csv", "a", lock=false) do f
                    writedlm(f, phi_new_out, ",")
                end
            else
                phi_t[:, :, t_index] = phi_new_out
            end
            if boundary == "neumann"
                mass_t[t_index] = calculate_mass(phi0, h2, nx / 2, ny / 2)
            else
                mass_t[t_index] = calculate_mass(phi0, h2, nx, ny)
            end
            E_t[t_index] = calculate_discrete_energy(phi_new_out, h2, epsilon2)
            D_t[t_index] = ch_r_error(r_new, phi_new, h2, C0, gamma0)

        end
        phi_prev = copy(phi_old)
        phi_old = copy(phi_new)
        r_old = copy(r_new)

    end
    # delta_mass_t = mass_t #.- mass_t[1]
    # E_t = E_t ./ E_t[1]
    return t_out, phi_t, mass_t, E_t
end