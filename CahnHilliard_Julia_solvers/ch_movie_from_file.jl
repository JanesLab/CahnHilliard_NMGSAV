using Plots, CSV, DataFrames

default(colorbar=false)  # Turn off colorbars globally for cleaner animation

function ch_movie_from_file(phi_file, t_out, ny; dtframes=1, filename="ch_movie", filetype="mp4", colorbar_type="default")
    """
    This function creates a red-white-blue video trajectory of chemical states using @animate.

    # Arguments
    - phi_file::String: File name for multidimensional array of chemical states.
    - t_out::Vector{Float64}: Vector of time steps corresponding to the third dimension of phi_t.
    - ny::Int: Number of mesh points in the y direction.

    # Optional Keyword Arguments
    - dtframes::Int: Number of dt time steps per frame. Default is 1. Note that this is a multiplier on dtout (i.e. if dtout = 10 and dtframes = 2, then every 20th timepoint is recorded in the movie). 
    - filename::String: Name of the movie file to be saved. Default is 'ch_movie'.
    - filetype::String: Format of the movie file ('mp4', 'gif'). Default is 'mp4'.
    - colorbar_type::String: Type of colorbar to be used ('default', 'initial_range'). Default is 'default'.

    # Output
    Saves a red-white-blue video.
    """

    # Validate file type
    valid_filetypes = ["mp4", "gif"]
    if !(filetype in valid_filetypes)
        error("Invalid filetype. Supported types are: mp4, gif")
    end

    # Initial color range
    initial_range = nothing

    anim = @animate for i in 1:dtframes:length(t_out)
        # Read the necessary frame from the file
        row_start = (i - 1) * ny + 1
        # row_end = i * ny
        df = DataFrame(CSV.File(phi_file; header=false, skipto=row_start, limit=ny))
        phi_temp = Matrix(df)

        # Determine color limits
        clims = if colorbar_type == "default"
            (-1, 1)
        elseif colorbar_type == "initial_range"
            initial_range === nothing && (initial_range = extrema(phi_temp))
            initial_range
        elseif colorbar_type == "variable"
            extrema(phi_temp)
        else
            error("Invalid colorbar_type: $colorbar_type")
        end

        # Create the heatmap
        title = @sprintf("%.3e", t_out[i])

        heatmap(phi_temp, y_flip=true, color=reverse(cgrad(:RdBu)), cbar=true, clims=clims, xlabel="", ylabel="", title="t = $(title)", xlim=(0, ny), ylim=(0, ny), aspect_ratio=:equal, dpi=300)
        if (i - 1) / (length(t_out) - 1) * 100 % 5 == 0
            println(@sprintf("%3.0f percent complete", (i - 1) / (length(t_out) - 1) * 100))
        end
    end #every dtframes <- redundant with above

    # Save animation
    if filetype == "mp4"
        mp4(anim, "$filename.mp4")
    elseif filetype == "gif"
        gif(anim, "$filename.gif")
    end

    println("Animation saved as $filename.$filetype")
end
