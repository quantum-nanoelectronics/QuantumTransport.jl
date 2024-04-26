using CSV
using DataFrames

# Adding Makie, etc. to the dependencies will not work on github
#using Makie #If this line is uncommented or if Makie is added to this package, github tests will fail

using CairoMakie



function generate_plot_makie(df, margin::Float64 = 0.0)
    # Check for invalid margin values
    if margin < 0
        throw(ArgumentError("Margin must be greater than 0"))
    end

    # Extract data directly from the DataFrame
    x = df[!, :X]
    y = df[!, :Y]
    z = df[!, :Z]
    d = df[!, :D]

    fig = Figure()
    ax = Axis3(
        fig[1, 1], 
        aspect = (1, 1, 1), 
        title = "Electron Density Plot",
        limits = (minimum(x)-margin, maximum(x)+margin, minimum(y)-margin, maximum(y)+margin, minimum(z)-margin, maximum(z)+margin),
        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z"
    )
    Colorbar(fig[1, 2], colormap = :viridis, flipaxis = false)
    
    scatter!(
        ax, 
        x, 
        y, 
        z, 
        color = d,
        markersize = d/5,
        colormap = :viridis
    )
    # display(fig)

    return fig
end



function plot_pos(df, nm::Float64)
    # Extract data directly from the DataFrame
    x = df[!, :X]
    y = df[!, :Y]
    z = df[!, :Z]
    d = df[!, :D]

    d = d / 100.0

    fig = meshscatter(x/nm, y/nm, z/nm, alpha=0.5,markersize=0.05)
    # fig = scatter(positions[1,:]/nm, positions[2,:]/nm, positions[3,:]/nm, type=:scatter,seriestype=:scatter, markersize = 4, c="black")
    # fig = plot(rand(5),rand(5))

    # display(fig)
    # use wait if display(fig) does not work
    # wait(display(fig))
    # readline()
    return fig
end