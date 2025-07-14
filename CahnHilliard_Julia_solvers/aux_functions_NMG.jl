using DataFrames
using DelimitedFiles
using LinearAlgebra
using BenchmarkTools
using StaticArrays
using Printf
import Random

#will need to ensure matrices are floats not ints
function convert_matrix(matrix)
    if eltype(matrix) == Float64
        return matrix  # Matrix is already of type Float64, no need to convert
    elseif eltype(matrix) == Int
        return convert(Matrix{Float64}, matrix)  # Convert to Matrix{Float64}
    else
        throw(ArgumentError("Unsupported matrix element type"))
    end
end

function calculate_mass(phi, h2, nx, ny)
    ave_mass = sum(phi) / (h2 * nx * ny)
    return ave_mass
end

function f(phi)
    fphi = (1 / 4) .* ((phi .^ 2) .- 1) .^ 2
    return fphi
end

function calculate_discrete_energy(phi, h2, epsilon2)
    a = h2 * sum(f(phi))
    sum_i = sum((phi[2:end, :] .- phi[1:end-1, :]) .^ 2)

    b = (epsilon2 / 2) * sum_i
    sum_j = sum((phi[:, 2:end] .- phi[:, 1:end-1]) .^ 2)

    c = (epsilon2 / 2) * sum_j
    E = a + b + c
    return E
end

function calculate_discrete_norm_energy(phi, phi0, h2, nx, ny, epsilon2)
    E0 = calculate_discrete_energy(phi0, h2, epsilon2)
    E = calculate_discrete_energy(phi, h2, epsilon2)
    return E / E0
end

function print_mat(file, matrix)
    open(file, "a", lock=false) do f
        for i = axes(matrix, 1)
            for j = axes(matrix, 2)
                Printf.@printf(f, "%16.15f ", matrix[i, j])
            end
            println(f)
        end
    end
end
