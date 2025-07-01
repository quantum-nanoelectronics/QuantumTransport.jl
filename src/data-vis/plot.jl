# New plotting file

# Adding Makie, etc. to the dependencies will not work on github
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail

# DPI, resolution?
const MAKIE_FIGURE_KWARGS = Dict(
    :size => (600, 400),
)

# font size, font, margin size?
const MAKIE_LINES_KWARGS = Dict([
    :linewidth => 1,
    :color => :black,
])

# colormap, color?
const MAKIE_AXIS_KWARGS = Dict([
    :xticksvisible => true, 
    :yticksvisible => true, 
    :xgridvisible => true, 
    :ygridvisible => true,
    :title => "(insert title here)",
])


# Function to call the appropriate function based on header value
function plot(readDir::String, filename::String ; kwargs...)
    println("Plotting data from ", filename)

    # Read the CSV file header
    df, metadata = get_data(readDir, filename)

    # metadata is in the format ["type of plot", n_entries, n_codomain (dimensionality of codomain), ...], where the ... are kwargs
    # f is the type_of_plot, stored as a symbol so that the function can be called directly
    type_of_plot, n_entries, dims_codomain, kwargs = metadata
    
    println("Type of plot: ", type_of_plot)
    println("Number of entries: ", n_entries)
    println("Dimensionality of codomain: ", dims_codomain)

    # Print key and value of additional keyword arguments
    println("All Keyword Arguments:")
    for (key, value) in kwargs
        println("  - $key: $value")
    end

    # Get the appropriate plotting function based on the type of plot
    f = getfield(@__MODULE__, type_of_plot)
    fig = f(df, dims_codomain; kwargs...)
    return fig
end

function ℝ_to_ℝ(df::DataFrame, n::Int; kwargs...)
    # TODO make this into a wrapper
    println("Plotting ℝ_to_ℝ data")
    filtered_figure_kwargs = merge(MAKIE_FIGURE_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_FIGURE_KWARGS) if haskey(kwargs, k)))
    println("Filtered Kwargs for Figure(): ", filtered_figure_kwargs)
    filtered_lines_kwargs = merge(MAKIE_LINES_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_LINES_KWARGS) if haskey(kwargs, k)))
    println("Filtered Kwargs for lines!(): ", filtered_lines_kwargs)
    filtered_axis_kwargs = merge(MAKIE_AXIS_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_AXIS_KWARGS) if haskey(kwargs, k)))
    println("Filtered Kwargs for Axis(): ", filtered_axis_kwargs)

    fig = Figure(; filtered_figure_kwargs...)
    xlabel = string(names(df)[1])
    ylabel = string(names(df)[2])

    if get(kwargs, :flip_axis, false)
        xlabel, ylabel = ylabel, xlabel
    end
    # methods available for Axis and lines! 
    # println(methods(lines!))
    # println(methods(Axis))

    ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, filtered_axis_kwargs...) 

    plot_ℝ_to_ℝ!(ax, df, n, filtered_lines_kwargs, kwargs) 

    # Display the figure
    #display(fig)

    filename = "R_to_R" * "_" * string(Dates.format(Dates.now(), "mm-dd_HH.MM.SS")) * ".png"
    output_file_path = string(joinpath(OUTPUT_DIR, filename))

    save(output_file_path, fig)

    return fig
end

function get_empty_plot(df::DataFrame ; kwargs...)
    # TODO make this into a wrapper
    println("Plotting ℝ_to_ℝ data")
    filtered_figure_kwargs = merge(MAKIE_FIGURE_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_FIGURE_KWARGS) if haskey(kwargs, k)))
    println("Filtered Kwargs for Figure(): ", filtered_figure_kwargs)
    filtered_lines_kwargs = merge(MAKIE_LINES_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_LINES_KWARGS) if haskey(kwargs, k)))
    println("Filtered Kwargs for lines!(): ", filtered_lines_kwargs)
    filtered_axis_kwargs = merge(MAKIE_AXIS_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_AXIS_KWARGS) if haskey(kwargs, k)))
    println("Filtered Kwargs for Axis(): ", filtered_axis_kwargs)

    fig = Figure(; filtered_figure_kwargs...)
    colors = distinguishable_colors(100)
    xlabel = string(names(df)[1])
    ylabel = string(names(df)[2])

    if get(kwargs, :flip_axis, false)
        xlabel, ylabel = ylabel, xlabel
    end

    axis = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, filtered_axis_kwargs...) 

    return fig, axis
end

function save_plot(fig::Figure, filename::String)
    # Save the figure
    output_file_path = string(joinpath(OUTPUT_DIR, filename))
    save(output_file_path, fig)
    println("Figure saved to ", output_file_path)
end

function add_to_ℝ_to_ℝ(ax, df::DataFrame, n::Int; kwargs...)
    # TODO make this into a wrapper
    println("Plotting ℝ_to_ℝ data")
    filtered_figure_kwargs = merge(MAKIE_FIGURE_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_FIGURE_KWARGS) if haskey(kwargs, k)))
    println("Filtered Kwargs for Figure(): ", filtered_figure_kwargs)
    filtered_lines_kwargs = merge(MAKIE_LINES_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_LINES_KWARGS) if haskey(kwargs, k)))
    println("Filtered Kwargs for lines!(): ", filtered_lines_kwargs)
    filtered_axis_kwargs = merge(MAKIE_AXIS_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_AXIS_KWARGS) if haskey(kwargs, k)))
    println("Filtered Kwargs for Axis(): ", filtered_axis_kwargs)

    plot_ℝ_to_ℝ!(ax, df, n, filtered_lines_kwargs, kwargs) 
end

function plot_ℝ_to_ℝ!(ax::Makie.Axis , df::DataFrame, n::Int, filtered_lines_kwargs::Dict, kwargs)
    domain = df[!, 1]
    colors = distinguishable_colors(100)
    xlabel = string(names(df)[1])
    ylabel = string(names(df)[2])

    if get(kwargs, :flip_axis, false)
        xlabel, ylabel = ylabel, xlabel
    end

    ax.xlabel = xlabel
    ax.ylabel = ylabel

    for i in 1:n
        codomain = df[!, i+1]  # get the (i+1)th column as the codomain
        # TODO not sure why we need a color_index here
        color_index = (i - 1) % length(colors) + 8  # cycle through colors list
        if get(kwargs, :flip_axis, true)
            domain, codomain = codomain, domain
        end
        lines!(ax, domain, codomain; color = colors[color_index], filtered_lines_kwargs...)
    end
end


# original function for reference
# function ℝ_to_ℝ(df::DataFrame, n::Int; kwargs...)
#     println("Plotting ℝ_to_ℝ data")
#     filtered_figure_kwargs = merge(MAKIE_FIGURE_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_FIGURE_KWARGS) if haskey(kwargs, k)))
#     println("Filtered Kwargs for Figure(): ", filtered_figure_kwargs)
#     filtered_lines_kwargs = merge(MAKIE_LINES_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_LINES_KWARGS) if haskey(kwargs, k)))
#     println("Filtered Kwargs for lines!(): ", filtered_lines_kwargs)
#     filtered_axis_kwargs = merge(MAKIE_AXIS_KWARGS, Dict(k => kwargs[k] for k in keys(MAKIE_AXIS_KWARGS) if haskey(kwargs, k)))
#     println("Filtered Kwargs for Axis(): ", filtered_axis_kwargs)

#     domain = df[!, 1]
#     fig = Figure(; filtered_figure_kwargs...)
#     colors = distinguishable_colors(100)
#     xlabel = string(names(df)[1])
#     ylabel = string(names(df)[2])

#     if get(kwargs, :flip_axis, false)
#         xlabel, ylabel = ylabel, xlabel
#     end
#     # methods available for Axis and lines! 
#     # println(methods(lines!))
#     # println(methods(Axis))

#     ax = Axis(fig[1, 1]; xlabel = xlabel, ylabel = ylabel, filtered_axis_kwargs...) 
    
#     # not sure why we need a color_index here
#     for i in 1:n
#         codomain = df[!, i+1]  # get the (i+1)th column as the codomain
#         color_index = (i - 1) % length(colors) + 8  # cycle through colors list
#         if get(kwargs, :flip_axis, true)
#             domain, codomain = codomain, domain
#         end
#         lines!(ax, domain, codomain; color = colors[color_index], filtered_lines_kwargs...)
#     end
#     # Display the figure
#     #display(fig)

#     filename = "R_to_R" * "_" * string(Dates.format(Dates.now(), "mm-dd_HH.MM.SS")) * ".png"
#     output_file_path = string(joinpath(OUTPUT_DIR, filename))

#     save(output_file_path, fig)

#     return fig
# end