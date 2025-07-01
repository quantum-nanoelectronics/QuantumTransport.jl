module TestScattering

using QuantumTransport
using Test
using Dates
using Colors
using CairoMakie

# Bypassing Driver for a sweep accross ϵ_rand_strengths

include(joinpath(INPUT_DIR, "AllInputs.jl"))
include("TestVisualizationFunctions.jl")


function get_latest_filename()
    # Get the latest filename from the directory
    return first(
        sort(
            filter(f -> startswith(f, "transmission"), readdir(OUTPUT_DIR));
            by = f -> stat(joinpath(OUTPUT_DIR, f)).mtime,
            rev = true
        )
    )
end


# you would just have to make the emptyy plot here
# then call the plot function to add to the plot

function fill_plot()
    figure, axis = nothing, nothing
    fig, ax = nothing, nothing
    
    for (i, ϵ) in enumerate(ϵ_values)
        p["ϵ_rand_strength"] = ϵ

        println("Running transport with random onsite disorder strength ϵ = $ϵ")

        # Run transport
        transport(p)

        # get the latest file
        filename = get_latest_filename()

        # plot(OUTPUT_DIR, filename)
        df, metadata = get_data(OUTPUT_DIR, filename)
        type_of_plot, n_entries, dims_codomain, kwargs = metadata

        if i == 1
            fig, ax = get_empty_plot(df ; kwargs...)
            figure, axis = fig, ax
        end

        kwargs = merge(kwargs, Dict(:color => colors[i], :linewidth => 1)) 

        # TODO make this a plot() call, need plot to accept kwargs and merge
        add_to_ℝ_to_ℝ(ax, df, dims_codomain; kwargs...)

    end
    return figure
end

p = runparams["transport"]
ϵ_values = [0.0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]

colors = distinguishable_colors(length(ϵ_values))
# colors = cgrad(:viridis, length(ϵ_values)).colors

fig = fill_plot()
filename = "scattering_plot" * "_" * string(Dates.format(Dates.now(), "mm-dd_HH.MM.SS")) * ".png"
save_plot(fig, filename)

end # module