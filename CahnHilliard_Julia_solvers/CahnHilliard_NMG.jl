# Julia implementation of a Multigrid solver for the Cahn-Hilliard equation with various optimizations for performance (particularly memory allocation)
# Author: Sarah Groves
# October 8, 2024

#=================================================
# Summary of Optimizations:
# [ ] Preallocate arrays and pass them to functions.
# [ ] Use in-place operations and @. for broadcasting.
# [ ] Add @inbounds and @simd for loop optimization.
# [ ] Use StaticArrays for small fixed-size arrays.
# [ ] Replace global constants with function arguments.
# [ ] Use @views to avoid array copying.
# [ ] Profile and benchmark to identify bottlenecks.
# [ ] Access arrays in a column-major order.
# [ ] Precompute and reuse frequently calculated constants.
# [ ] Consider SparseArrays for sparse matrices.
# [ ] Ensure type stability with type annotations.
=================================================#

using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random
include("nmg_solver.jl")
include("ch_initialization.jl")
include("aux_functions_NMG.jl")

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

    OUTPUT
        t_out = Time corresponding to the dt time step outputs.
        phi_t = Multidimensional array of phi over t_out.
        delta_mass_t = Vector of mass change over t_out.
        E_t = Vector of relative energy over t_out.
"""
function CahnHilliard_NMG(phi0; t_iter=1e3, dt=2.5e-5, solver_iter=1e4, tol=1e-5, dt_out=10, m=8, epsilon2=NaN, boundary="periodic", c_relax=2, domain=[1 0 1 0], printres=false, printphi=false, pathname="cd")
    nx, ny = size(phi0)
    xright, xleft, yright, yleft = domain
    # ny = nx
    n_level::Int = trunc(log(nx) / log(2.0) + 0.1)  # original c code uses natural log too
    h2 = ((xright - xleft) / nx) * ((yright - yleft) / ny) #Define mesh size
    if isnan(epsilon2)
        epsilon2 = h2 * m^2 / (2 * sqrt(2) * atanh(0.9))^2 #Define ϵ^2 if not prespecified
    else
        m = sqrt((epsilon2 * (2 * sqrt(2) * atanh(0.9))^2) / h2) #Else overwrite m
    end

    mu = zeros(Float64, nx, ny)
    phi_new = copy(phi0)
    phi_old = copy(phi0)

    downsampled = nx * ny * t_iter / dt_out > 1e9 #Logical index for the need to downsample
    n_timesteps = floor(Int64, t_iter / dt_out)

    # Save or print initial phi
    if printphi
        mass_t = zeros(Float64, n_timesteps + 1, 1)
        E_t = zeros(Float64, n_timesteps + 1, 1)
        t_out = 0:dt_out*dt:(n_timesteps)*dt*dt_out
        if pathname == "cd"
            pathname = pwd()
        end
        open("$(pathname)phi.csv", "w", lock=false) do f
            writedlm(f, phi0, ",")
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
        t_out = 0:dt_out*dt:(n_timesteps)*dt*dt_out
        phi_t = zeros(Float64, nx, ny, n_timesteps + 1)
        phi_t[:, :, 1] = phi0
    end

    mass_t[1] = calculate_mass(phi0, h2, nx, ny)
    E_t[1] = calculate_discrete_energy(phi0, h2, epsilon2)
    if printres
        println("Saving squared residuals per iteration to file in the output directory\n")
    end
    rr = zeros(Float64, nx, ny)
    for it in 1:t_iter
        nmg_solver!(rr, phi_old, phi_new, mu, nx, ny, dt, solver_iter, tol, c_relax, xright, xleft, yright, yleft, epsilon2, n_level, boundary, printres=printres)
        phi_old = copy(phi_new)
        if it / t_iter * 100 % 5 == 0
            println(@sprintf("%3.0f percent complete", it / t_iter * 100))
        end
        if it % dt_out == 0
            t_index = floor(Int64, it / dt_out) + 1
            if printphi
                open("$(pathname)phi.csv", "a", lock=false) do f
                    writedlm(f, phi_new, ",")
                end
            else
                phi_t[:, :, t_index] = phi_new
            end
            mass_t[t_index] = calculate_mass(phi_new, h2, nx, ny)
            E_t[t_index] = calculate_discrete_energy(phi_new, h2, epsilon2)
        end
    end
    delta_mass_t = mass_t #.- mass_t[1]
    # E_t = E_t ./ E_t[1]
    return t_out, phi_t, delta_mass_t, E_t
end