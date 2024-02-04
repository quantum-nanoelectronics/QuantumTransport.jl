using QuantumTransport
using Test

# test the files in the io folder

const SAMPLE_POSITIONS = [
    [1.0, 0.0, 0.0],
    [2.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
    [4.0, 0.0, 0.0],
    [5.0, 0.0, 0.0],
    [6.0, 0.0, 0.0],
    [7.0, 0.0, 0.0],
    [8.0, 0.0, 0.0],
    [9.0, 0.0, 0.0],
    [10.0, 0.0, 0.0]
]

baseDir = abspath(joinpath(@__DIR__, ".."))
ioDir = joinpath(baseDir, "data-output")
df, meta = generate_csv(SAMPLE_POSITIONS)

@test save_csv(ioDir, df, meta)

println("-Writing-")
println("DataFrame: ")
println(df)
println("Metadata: ")
println(meta)

vals = get_data(ioDir)
@test !isnothing(vals[1])

println("-Reading-")
println("DataFrame: ")
println(vals[1])
println("Metadata: ")
println(vals[2])
