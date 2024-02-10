# test the files in the data-vis folder

using QuantumTransport
using Test

# contains constants used in make-cnt-data and other data visualization
const nm = 1E-9
const δ₀ = 0.142*nm
const a = 0.246*nm

margin = 0.0

# Assuming the required file exists in the io folder
baseDir = abspath(joinpath(@__DIR__, ".."))
ioDir = joinpath(baseDir, "data-output")


vals = get_data(ioDir)

df = vals[1]
metaData = vals[2]
@test !isnothing(df)

println("-Reading-")
println("DataFrame: ")
println(first(df, 5))
println("Metadata: ")
println(metaData)

@test plot_pos(df, nm)
@test generate_plot_makie(df, margin)
