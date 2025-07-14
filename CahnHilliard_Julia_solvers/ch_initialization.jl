using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random

# function initialize_geometric_CPC(nx, ny; CPC_width=20, cohesin_width=4)
#     phi = zeros(Float64, nx, ny)
#     # CPC_width = 20
#     # cohesin_width = 4
#     @simd for i in 1:nx
#         for j in 1:ny
#             if i > round(nx / 2) - CPC_width && i < round(nx / 2) + CPC_width
#                 if j > round(ny / 2) - CPC_width && j < round(ny / 2) + CPC_width
#                     phi[i, j] = 1
#                 elseif i > round(nx / 2) - cohesin_width && i < round(nx / 2) + cohesin_width
#                     phi[i, j] = 1
#                 else
#                     phi[i, j] = -1
#                 end
#             else
#                 phi[i, j] = -1
#             end
#         end
#     end
#     return phi
# end

function initialize_round_CPC_um(nx, ny; CPC_radius=0.173, cohesin_width=0.2, domain_width=3.2, c1=-1.0, c2=1.0)
    CPC_radius_grid_points = nx * CPC_radius / domain_width
    cohesin_width_grid_points = nx * cohesin_width / domain_width
    cohesin_half_width = cohesin_width_grid_points / 2
    # Create an empty matrix filled with -1
    phi = fill(c1, nx, ny)

    # Define the center of the matrix
    center = (nx) / 2

    # Loop through each element of the matrix
    for i in 1:nx
        for j in 1:ny
            # Calculate the distance from the center
            distance = norm([i - center, j - center])

            # Check if the distance is less than or equal to CPC_width
            if distance <= CPC_radius_grid_points
                phi[i, j] = c2
            elseif i > round((nx) / 2) - cohesin_half_width && i < round((nx) / 2) + cohesin_half_width
                phi[i, j] = c2
            end
        end
    end
    return phi
end

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

function initialization_from_function(nx, ny; R0=0.1, epsilon=0.01)
    phi = zeros(Float64, nx, ny)
    h = 1 / nx
    x = h .* (0:nx-1)
    y = h .* (0:ny-1)
    xx, yy = meshgrid(x, y)
    R = @.sqrt((xx - 0.5)^2 + (yy - 0.5)^2)
    # eps_c = epsilon
    # delta = eps_c * sqrt(2)
    # psi0 = 0.5 * (1 .+ @.tanh((R0 .- R) / (2 * delta)))
    # phi = 2 .* psi0 .- 1    # psi0=(phi0+1)/2
    phi = @.tanh((R0 .- R) / (sqrt(2) * epsilon))
    return phi
end

function initialization_random(nx, ny)
    return (1 .- 2 .* rand(nx, ny))
end

function initialization_spinodal(nx, ny)
    return (rand([-1.0, 1.0], nx, ny))
end

function initialization_from_file(file, nx, ny; delim=',', transpose_matrix=false)
    phi = readdlm(file, delim, Float64)
    if size(phi) != (nx, ny)
        print("Warning: phi from file is wrong size: $(size(phi)) Expected: $(nx), $(ny)")
    end
    if transpose_matrix
        phi = transpose(phi)
    end
    return phi
end

function ch_initialization(nx, ny; method="spinodal", initial_file="", delim=",", R0=0.1, epsilon=0.01, cohesin_width=0.1, CPC_radius=0.2, domain_width=3.2, c1=-1.0, c2=1.0)
    if method == "random"
        phi0 = initialization_random(nx, ny)
    elseif method == "droplet"
        phi0 = initialization_from_function(nx, ny, R0=R0, epsilon=epsilon)
    elseif method == "geometric"
        phi0 = initialize_round_CPC_um(nx, ny, CPC_radius=CPC_radius, cohesin_width=cohesin_width, domain_width=3.2, c1=-1.0, c2=1.0)
    elseif method == "file"
        phi0 = initialization_from_file(initial_file, nx, ny, delim=delim)
    elseif method == "spinodal"
        phi0 = initialization_spinodal(nx, ny)
    else
        println("Warning: initialize must be one of [random, droplet, geometric, file, spinodal].")
    end
    return phi0
end
