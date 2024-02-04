using QuantumTransport
using Test

# test the files in the io folder

baseDir = abspath(joinpath(@__DIR__, ".."))
ioDir = joinpath(baseDir, "data-output")
@test generate_csv(ioDir)
@test !isnothing(get_data(ioDir))

