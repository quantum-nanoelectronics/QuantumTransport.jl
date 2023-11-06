using CSV
using DataFrames
using GLMakie

#include("../io/read-positions.jl")
include("constants.jl")



function generate_plot_makie(margin::Float64 = 0.0)
    # Check for invalid margin values
    if margin < 0
        throw(ArgumentError("Margin must be greater than 0"))
    end

    x, y, z, d = get_data()
    zs = LinRange(0, 3, length(x))


    fig = Figure()
    ax = Axis3(
        fig[1, 1], 
        aspect = (1, 1, 1), 
        title = "Electron Density Plot",
        limits = (minimum(x)-margin, maximum(x)+margin, minimum(y)-margin, maximum(y)+margin, minimum(z)-margin, maximum(z)+margin),
        xlabel = "X",
        ylabel = "Y",
        zlabel = "Z",

        )
    Colorbar(fig[1, 2], colormap = :viridis, flipaxis = false) 
    
    GLMakie.scatter!(
        ax, 
        x, 
        y, 
        z, 
        color = d,
        markersize = d/5,
        colormap = :viridis, 
        )
    display(fig)

    return true
end


function plot_pos()
    x, y, z, d = get_data()
    d = d / 100.0

    fig = meshscatter(x/nm, y/nm, z/nm, alpha=0.5,markersize=0.05)
    # fig = scatter(positions[1,:]/nm, positions[2,:]/nm, positions[3,:]/nm, type=:scatter,seriestype=:scatter, markersize = 4, c="black")
    # fig = plot(rand(5),rand(5))

    display(fig)
    # use wait if display(fig) does not work
    # wait(display(fig))
    # readline()
    return true
end

# generate_plot_makie()
# plot_pos()























###
### Archived code below
###


# using CairoMakie

# old function with Plots library
# function generate_plot()
#     x, y, z, d = get_data()
#     d = d / 100.0


#     colors = ifelse.(d .< 50, :orange, :blue)
#     color_scale = cgrad([RGB(0.8, 0.8, 1), RGB(0, 0, 0.5)], d)


#     plt3d= Plots.plot(x, y, z,seriestype=:scatter, color = color_scale, marker_z = d, markersize = 3, title="Electron Positions", xlabel="X", ylabel="Y", zlabel="Z")
#     #plt3d = Plots.scatter(x, y, z, color = d, marker_z = d, colorbar = (title = "Color"), c = color_scale)
    
#     display(plt3d)
# end

# function generate_plot_makie(margin)

#     # x = cos.(1:0.5:20)
#     # y = sin.(1:0.5:20)
#     # z = LinRange(0, 3, length(x))

#     #  s = CairoMakie.meshscatter(
#     #     x, 
#     #     y, 
#     #     z, 
#     #     markersize = 0.1,
#     #     color = d,

#     # )
#     # s
#     #limits = (minimum(x)-margin, maximum(x)+margin, minimum(y)-margin, maximum(y)+margin, minimum(z)-margin, maximum(z)+margin))

# end

#generate_plot()

# generate_plot_makie(5)


# y = rand(100)
# plot(seriestype=:scattery, zcolor = abs.(y .- 0.5), m = (:heat, 0.8, Plots.stroke(1, :green)), ms = 10 * abs.(y .- 0.5) .+ 4, lab = "grad")
# scatter!(y, zcolor = abs.(y .- 0.5), m = (:heat, 0.8, Plots.stroke(1, :green)), ms = 10 * abs.(y .- 0.5) .+ 4, lab = "grad")



# default(legend = false)
# x = y = range(-5, 5, length = 40)
# zs = zeros(0, 40)
# n = 100

# @gif for i in range(0, stop = 2π, length = n)
#     f(x, y) = sin(x + 10sin(i)) + cos(y)

#     # create a plot with 3 subplots and a custom layout
#     l = @layout [a{0.7w} b; c{0.2h}]
#     p = plot(x, y, f, st = [:surface, :contourf], layout = l)

#     # induce a slight oscillating camera angle sweep, in degrees (azimuth, altitude)
#     plot!(p[1], camera = (10 * (1 + cos(i)), 40))

#     # add a tracking line
#     fixed_x = zeros(40)
#     z = map(f, fixed_x, y)
#     plot!(p[1], fixed_x, y, z, line = (:black, 5, 0.2))
#     vline!(p[2], [0], line = (:black, 5))

#     # add to and show the tracked values over time
#     global zs = vcat(zs, z')
#     plot!(p[3], zs, alpha = 0.2, palette = cgrad(:blues).colors)
# end





