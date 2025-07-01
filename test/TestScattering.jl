module TestScattering

using QuantumTransport
using Test
using Dates
using Colors

# Bypassing Driver for a sweep across ϵ_rand_strengths

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

function fill_plot()
    title = "Scattering Transport with Random Onsite Disorder"
    get_empty_plot(; size = (600, 400), title = title, xlabel = "E (eV)", ylabel = "T (e²/h)")
    fig, ax = get_empty_plot(; )
    
    for (i, ϵ) in enumerate(ϵ_values)
        p["ϵ_rand_strength"] = ϵ

        println("Running transport with random onsite disorder strength ϵ = $ϵ")

        # Run transport
        transport(p)

        # get the latest file
        filename = get_latest_filename()

        # Plot the transmission data
        plot(OUTPUT_DIR, filename; axis = ax, color = colors[i], linewidth = 1, scattering = p["scattering"], title = title)

    end
    return fig
end

p = runparams["transport"]
ϵ_values = [0.0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0]

colors = distinguishable_colors(length(ϵ_values))
# colors = cgrad(:viridis, length(ϵ_values)).colors

fig = fill_plot()
filename = "scattering_plot" * "_" * string(Dates.format(Dates.now(), "mm-dd_HH.MM.SS")) * ".png"
save_plot(fig, filename)

end # module