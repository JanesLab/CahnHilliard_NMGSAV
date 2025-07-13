using DelimitedFiles
using Dates
using BenchmarkTools
using Printf
include("CahnHilliard_NMG.jl")
include("CahnHilliard_SAV.jl")
include("ch_movie_from_file.jl")

indir = "./IC/"
outdir = "./output/output_julia/"
n_relax = 4
GridSize = 128
total_time = 0.02
h = 1 / GridSize
m = 8
epsilon = m * h / (2 * sqrt(2) * atanh(0.9))

boundary = "periodic"
printphi = true
phi0 = initialization_from_file("$(indir)initial_phi_$(GridSize)_smooth_n_relax_$n_relax.csv", GridSize, GridSize, delim=',', transpose_matrix=false)
dt_out = 10

# for dt in [5.5e-6, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4]
for GridSize = [64, 128, 256]
    dt = 5.5e-6
    if dt == 5.5e-6
        max_it = 2000
    else
        max_it = floor(Int64, total_time / dt)
    end

    ################
    # NMG
    ################
    method = "NMG"
    pathname = @sprintf("%s%s_Julia_%d_dt_%.2e_Nx_%d_n_relax_%d_", outdir, method, max_it, dt, GridSize, n_relax)
    date_time = now()
    result, elapsed_time, mem_allocated, gc_time, memory_counters = @timed CahnHilliard_NMG(phi0; t_iter=max_it, dt=dt, dt_out=dt_out, m=m, boundary=boundary, printphi=printphi, pathname=pathname)
    open("./Job_specs.csv", "a", lock=false) do f
        writedlm(f, [Dates.format(date_time, "mm/dd/yyyy HH:MM:SS") "spinodal_decomp_smoothed_dt_scan" "Julia" method GridSize epsilon dt 1e-5 max_it 1e4 pathname elapsed_time mem_allocated / 1e6], ",")
    end
    t_out = result[1]
    phi_t = result[2]
    mass_t = result[3]
    E_t = result[4]
    open("$(pathname)t_out.csv", "w", lock=false) do f
        writedlm(f, t_out, " ")
    end
    open("$(pathname)mass.csv", "w", lock=false) do f
        writedlm(f, mass_t, " ")
    end
    open("$(pathname)energy.csv", "w", lock=false) do f
        writedlm(f, E_t, " ")
    end
    # t_out = readdlm("$(pathname)t_out.csv", ',', Float64)
    # ch_movie_from_file("$(pathname)phi.csv", t_out, 128; dtframes=1, filename="$(pathname)movie", filetype="mp4", colorbar_type="default")

    ########
    # SAV
    ########
    method = "SAV"
    pathname = @sprintf("%s%s_SAV_%d_dt_%.2e_Nx_%d_n_relax_%d_", outdir, method, max_it, dt, GridSize, n_relax)

    date_time = now()
    result, elapsed_time, mem_allocated, gc_time, memory_counters = @timed CahnHilliard_SAV(phi0; t_iter=max_it, dt=dt, dt_out=dt_out, m=m, boundary=boundary, printphi=printphi, pathname=pathname)
    open("./Job_specs.csv", "a", lock=false) do f
        writedlm(f, [Dates.format(date_time, "mm/dd/yyyy HH:MM:SS") "spinodal_decomp_smoothed_dt_scan" "Julia" method GridSize epsilon dt 1e-5 max_it 1e4 pathname elapsed_time mem_allocated / 1e6], ",")
    end

    t_out = result[1]
    phi_t = result[2]
    mass_t = result[3]
    E_t = result[4]
    open("$(pathname)t_out.csv", "w", lock=false) do f
        writedlm(f, t_out, " ")
    end
    open("$(pathname)mass.csv", "w", lock=false) do f
        writedlm(f, mass_t, " ")
    end
    open("$(pathname)energy.csv", "w", lock=false) do f
        writedlm(f, E_t, " ")
    end
    # t_out = readdlm("$(pathname)t_out.csv", ',', Float64)
    # ch_movie_from_file("$(pathname)phi.csv", t_out, 128; dtframes=1, filename="$(pathname)movie", filetype="mp4", colorbar_type="default")
end