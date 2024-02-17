# test the files in the data-vis folder

using QuantumTransport
using Test
using CairoMakie

# contains constants used in make-cnt-data and other data visualization
const nm = 1E-9
const δ₀ = 0.142*nm
const a = 0.246*nm

margin = 0.0
filename = "scatterplot.csv"

# Assuming the required file exists in the io folder
baseDir = abspath(joinpath(@__DIR__, ".."))
ioDir = joinpath(baseDir, "data-output")


vals = get_data(ioDir, filename)

df = vals[1]
metaData = vals[2]
# @test !isnothing(df)

println("-Reading-")
println("DataFrame: ")
println(first(df, 5))
println("Metadata: ")
println(metaData)


# Obtain the figure object
figure = generate_plot_makie(df, margin)
@test figure isa Figure
# Save the figure
file_path = joinpath(ioDir, "figure_from_generate_plot_makie.png")
save(file_path, figure)
@test isfile(file_path)

# Obtain the figure object
figure = plot_pos(df, nm)
@test !isnothing(figure)
# Save the figure
file_path = joinpath(ioDir, "figure_from_plot_pos.png")
save(file_path, figure)
@test isfile(file_path)
