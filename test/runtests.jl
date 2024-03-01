using QuantumTransport
using Test

"""
This is the main test file that is called first when running any tests.
All other dependencies for tests besides the QuantumTransport package and the Test package should be included in this file.
"""


@testset "QuantumTransport.jl Tests" begin

    # This will include and run the tests from `hello_world.jl`
    println("---------------RUNNING SAMPLE TESTS---------------")
    @testset "Sample Test" begin
        include("hello-world.jl")
        include("column-major.jl")
    end

    println("----------RUNNING MATRIX INVERSION TESTS----------")
    @testset "Matrix Inversion Test" begin
        using LinearAlgebra
        using SparseArrays
        include("matrix.jl")
    end

    println("------------RUNNING INPUT OUTPUT TESTS------------")
    @testset "Input Output Test" begin
        using DataFrames
        using Random
        include("io.jl")
    end

    println("------------RUNNING SELF ENERGIES TESTS------------")
    @testset "Self Energies Test" begin
        include("self-energies.jl")
    end

    # cannot fully test Data Visualization because GitHub does not have a GPU
    # non interactive image plots generated, however we want to generate interactive plots
    @testset "Data Visualization Test" begin
        using CairoMakie
        include("visualization.jl")
    end

end
