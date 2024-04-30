#New plotting file

# Adding Makie, etc. to the dependencies will not work on github
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail

# Function to call the appropriate function based on header value
function call_function_based_on_header(readDir::String, filename::String)
    # Read the CSV file header
    df, metadata = get_data(readDir, filename)

    # Get the header of the DataFrame
    if metadata[1] == "ℝ→ℝ"
        fig = ℝtoℝ(df, parse(Int, metadata[3]), parse(Bool, metadata[5]), String(metadata[4]::SubString))
    elseif metadata[1] == "ℝ³→ℝ"
        fig = ℝ³toℝ(df, parse(Int, metadata[3]), parse(Bool, metadata[5]), String(metadata[4]::SubString))
    else
        print("no match")
    end

    return fig
end


function ℝtoℝ(df::DataFrame, n::Int, axisflag::Bool, title::String, plot_size=(600, 400), line_scale=2)
    # Read the CSV file
      

    # Extract the domain and codomains
    domain = df[!, 1]  # first column for the domain R

    # Create the plot with the specified size
    fig = Figure(size = plot_size)
    colors = distinguishable_colors(100)
    if axisflag
        ax = Axis(fig[1, 1], title = title, xlabel = string(names(df)[2]), ylabel = string(names(df)[1])) 
    else
        ax = Axis(fig[1, 1], title = title, xlabel = string(names(df)[1]), ylabel = string(names(df)[2])) 
    end
     
    for i in 1:n
        codomain = df[!, i+1]  # get the (i+1)th column as the codomain
        color_index = (i - 1) % length(colors) + 8  # cycle through colors list
        if axisflag
            lines!(ax, codomain, domain, color = colors[color_index], linewidth = 1.5 * line_scale)
        else
            lines!(ax, domain, codomain, color = colors[color_index], linewidth = 1.5 * line_scale)
        end
    end

    # Display the figure
    #display(fig)
    output_file_path = joinpath(OUTPUT_DIR, "R_to_R.png")
    save(output_file_path, fig)

    return fig
end

function R_to_C(csv_file, plot_size=(600, 400), line_scale=1)
    
    df = CSV.File(csv_file) |> DataFrame

    # Function to map phase to inferno colormap
    function phase_to_inferno(phase)
        # Normalize phase to the range [0, 1]
        normalized_phase = (phase % (2pi)) / (2pi)
        # Map normalized phase to inferno colormap
        color_index = max(1, ceil(Int, normalized_phase * length(ColorSchemes.inferno)))
        return ColorSchemes.inferno[color_index]
    end

    # Extract the domain and the real and imaginary parts of the codomain
    domain = df[!, 1]  # First column for the domain R
    complex_num = parse.(Complex{Float64}, df[!, 2])  # Assuming second column has complex numbers in string format
    real_part = real.(complex_num)
    imag_part = imag.(complex_num)

    # Compute magnitude and phase of the complex numbers
    magnitude = sqrt.(real_part.^2 + imag_part.^2)
    phase = atan.(imag_part, real_part)  # atan2 is used to compute the angle

    # Normalize magnitudes for line widths
    max_magnitude = maximum(magnitude)
    normalized_magnitude = magnitude ./ max_magnitude * 5 * line_scale  # Adjust and scale line width

    # Create the plot with the specified size
    fig = Figure(size = plot_size)
    ax = fig[1, 1] = Axis(fig, title="Magnitude and Phase Plot", xlabel="R (domain)", ylabel="|C| (magnitude)")

    # Add the line segments with color based on phase and width based on magnitude
    for i in 1:length(domain) - 1
        # Determine the color based on the phase
        color = phase_to_inferno(phase[i])
        # Plot a segment of the line with the determined color and scaled line width
        lines!(ax, domain[i:i+1], magnitude[i:i+1], color=color, linewidth=normalized_magnitude[i])
    end

    # Display the figure
    #display(fig)
    output_file_path = joinpath(output_dir, "R_to_C.png")
    save(output_file_path, fig)

    return fig
end

function R3_to_Rplus(csv_file, plot_size=(600, 400), marker_scale=1)
      

    # Extract the coordinates and complex number values
    x = df[!, 1]  # First column for the x coordinate
    y = df[!, 2]  # Second column for the y coordinate
    z = df[!, 3]  # Third column for the z coordinate
    real_part = df[!, 4]


    # Determine marker sizes based on the magnitude of the real part, scaled by the marker_scale parameter
    marker_sizes = abs.(real_part) * marker_scale
    max_size = maximum(marker_sizes)
    normalized_sizes = marker_sizes ./ max_size   # Apply the scale factor for sizes

    # Create the 3D scatter plot with the specified size
    fig = Figure(size = plot_size)
    ax = fig[1, 1] = Axis3(fig, title=title, xlabel=string(names(df)[1]), ylabel=string(names(df)[2]), zlabel=string(names(df)[3]))

    # Normalize real part values for color mappings
    normalized_values = (real_part .- minimum(real_part)) ./ (maximum(real_part) - minimum(real_part))

    # Add the scatter points to the plot, using the seismic colormap for colors based on the normalized real part values
    scatter!(ax, x, y, z, markersize=normalized_sizes, color=normalized_values, colormap=:viridis, marker_z=normalized_values, strokewidth = 0)
    output_file_path = joinpath(OUTPUT_DIR, "R3_to_R.png")
    save(output_file_path, fig)

    return fig
end

function ℝ³toℝ(df::DataFrame, n::Int, axisflag::Bool, title::String, plot_size=(1400, 800), marker_scale=20)
    # Read the CSV file
      

    # Extract the coordinates and complex number values
    x = df[!, 1]  # First column for the x coordinate
    y = df[!, 2]  # Second column for the y coordinate
    z = df[!, 3]  # Third column for the z coordinate
    real_part = df[!, 4]


    # Determine marker sizes based on the magnitude of the real part, scaled by the marker_scale parameter
    marker_sizes = abs.(real_part) 
    max_size = maximum(marker_sizes)
    normalized_sizes = marker_sizes ./ max_size .* marker_scale # Apply the scale factor for sizes

    # Create the 3D scatter plot with the specified size
    fig = Figure(size = plot_size)
    ax = fig[1, 1] = Axis3(fig, title=title, xlabel=string(names(df)[1]), ylabel=string(names(df)[2]), zlabel=string(names(df)[3]))


    # Add the scatter points to the plot, using the seismic colormap for colors based on the normalized real part values
    scatter_plot = scatter!(ax, x, y, z, markersize=normalized_sizes, color=normalized_sizes, colormap=:viridis, strokewidth = 0)
    cb = Colorbar(fig[1, 2], limits = (0, maximum(marker_sizes)), label=string(names(df)[4]), )
    output_file_path = joinpath(OUTPUT_DIR, "R3_to_R.png")
    save(output_file_path, fig)
    
    return fig
end

function R3_to_C(csv_file, plot_size=(600, 400), marker_scale=1)
    
    df = CSV.File(csv_file) |> DataFrame

    # Extract the coordinates and complex number components
    x = df[!, 1]  # X coordinate
    y = df[!, 2]  # Y coordinate
    z = df[!, 3]  # Z coordinate
    complex_num = parse.(Complex{Float64}, df[!, 4])
    real_part = real.(complex_num)  # Real part of the complex number
    imag_part = imag.(complex_num)  # Imaginary part of the complex number

    # Calculate the phase of the complex numbers
    phases = atan.(imag_part, real_part)  # atan2 for phase calculation

    # Normalize phase to the range [0, 1]
    normalized_phases = (phases .% (2pi)) / (2pi)

    # Normalize sizes for better visualization and apply marker scaling
    magnitudes_squared = real_part.^2 + imag_part.^2
    max_magnitude = maximum(magnitudes_squared)
    normalized_sizes = (sqrt.(magnitudes_squared ./ max_magnitude) .* 30) * marker_scale  # Apply sqrt, then scale

    # Create the plot with the specified size
    fig = Figure(size = plot_size)
    ax = fig[1, 1] = Axis3(fig, title="3D Complex Function Visualization", xlabel="X", ylabel="Y", zlabel="Z")

    # Add the scatter points to the plot, using the inferno colormap for colors based on the phase
    scatter!(ax, x, y, z, color=normalized_phases, colormap=:inferno, markersize=normalized_sizes)

    # Display the figure
    #display(fig)
    output_file_path = joinpath(output_dir, "R3_to_C.png")
    save(output_file_path, fig)

    return fig
end

function codomain_to_HSL(codomain_y1, codomain_y2, codomain_y3, R_max)
    # Calculate Hue
    H = (atan(codomain_y1, codomain_y2) + π) * (180 / π)  # Correct the angle and convert radians to degrees
    H = mod(H, 360)  # Ensure H is within the [0, 360] range
    
    # Calculate Saturation
    S = sqrt(codomain_y1^2 + codomain_y2^2 + codomain_y3^2) / R_max
    S = clamp(S, 0, 1)  # Ensure S is within the [0, 1] range
    
    # Calculate Lightness
    L = codomain_y3 / R_max
    L = clamp(L, 0, 1)  # Ensure L is within the [0, 1] range
    
    return HSL(H, S, L)
end

function R3_to_R3(csv_file, plot_size=(600, 400), marker_scale=1)
    
    df = CSV.File(csv_file) |> DataFrame

    # Extract the domain and codomain
    domain_x, domain_y, domain_z = df[!, 1], df[!, 2], df[!, 3]
    codomain_x, codomain_y, codomain_z = df[!, 4], df[!, 5], df[!, 6]

    # Normalize or scale codomain_x for hue (0 to 360), codomain_y and codomain_z for saturation and value (0 to 1)

    # Calculate the marker sizes based on the magnitude of the spin density vectors and scale them
    marker_sizes = sqrt.(codomain_x.^2 + codomain_y.^2 + codomain_z.^2)
    max_marker_size = maximum(marker_sizes)
    normalized_marker_sizes = (sqrt.(marker_sizes ./ max_marker_size) .* 30) * marker_scale

    # Create the plot with the specified size
    fig = Figure(size = plot_size)
    ax = fig[1, 1] = Axis3(fig, title="3D Spin Density Visualization", xlabel="X", ylabel="Y", zlabel="Z")

    # Add the scatter points to the plot, with size and color
    colors = [codomain_to_HSL(codomain_x[i], codomain_y[i], codomain_z[i], max_marker_size) for i in 1:length(codomain_x)]
    scatter!(ax, domain_x, domain_y, domain_z, color=colors, markersize=normalized_marker_sizes)

    # Display the figure
    #display(fig)
    output_file_path = joinpath(output_dir, "R3_to_R3.png")
    save(output_file_path, fig)

    return fig
end

function R3(csv_file, hex_color, plot_size=(600, 400), marker_scale=1)
    # Read the CSV file
    df = CSV.File(csv_file) |> DataFrame

    # Extract the coordinates
    x = df[!, 1]  # X coordinate
    y = df[!, 2]  # Y coordinate
    z = df[!, 3]  # Z coordinate

    # Create the plot with the specified size
    fig = Figure(size = plot_size)
    ax = fig[1, 1] = Axis3(fig, title="3D Scatter Plot with Uniform Color", xlabel="X", ylabel="Y", zlabel="Z")

    # Add the scatter points to the plot with the specified hex color and scaled marker size
    scatter!(ax, x, y, z, color=hex_color, markersize=marker_scale)
    #display(fig)
    # Display the figure
    output_file_path = joinpath(output_dir, "R3.png")
    save(output_file_path, fig)
end

function R2_to_Rplus(csv_file, plot_size=(600, 400))
    
    # Read the CSV file into a DataFrame
    df = CSV.File(csv_file) |> DataFrame

    # Extract x, y coordinates, and values
    x_coords = unique(df[!, 1])
    y_coords = unique(df[!, 2])
    values_matrix = reshape(df[!, 3], length(y_coords), length(x_coords))

    # Create the plot with the specified size
    fig = Figure(size = plot_size)
    ax = fig[1, 1] = Axis(fig, title="2D Heatmap using Inferno Colormap")

    # Add the heatmap to the plot
    heatmap!(ax, x_coords, y_coords, values_matrix, colormap=:inferno)

    # Display the figure
    #display(fig)
    output_file_path = joinpath(output_dir, "R2_to_Rplus.png")
    save(output_file_path, fig)

    return fig
end

function R2_to_R3(csv_file, plot_size=(600, 400), marker_scale=1)
    
    # Read the CSV file into a DataFrame
    df = CSV.File(csv_file) |> DataFrame

    # Extract the coordinates
    x = df[!, 1]  # X coordinate
    y = df[!, 2]  # Y coordinate

    # Extract the RGB color values
    codomain_x, codomain_y, codomain_z = df[!, 3], df[!, 4], df[!, 5]

    # Normalize the RGB values using the range method
    marker_sizes = sqrt.(codomain_x.^2 + codomain_y.^2 + codomain_z.^2)
    max_marker_size = maximum(marker_sizes)

    # Create the plot with the specified size
    fig = Figure(size = plot_size)
    ax = fig[1, 1] = Axis(fig, title="2D Scatter Plot with Normalized HSL", xlabel="X", ylabel="Y")

    # Combine the normalized RGB values into a single color per point
    colors = [codomain_to_HSL(codomain_x[i], codomain_y[i], codomain_z[i], max_marker_size) for i in 1:length(codomain_x)]

    # Add a scatter plot to the plot, scaling the marker size
    scatter!(ax, x, y, color=colors, marker=:circle, markersize=20 * marker_scale)

    # Display the figure
    #display(fig)
    output_file_path = joinpath(output_dir, "R2_to_R3.png")
    save(output_file_path, fig)

    return fig
end

