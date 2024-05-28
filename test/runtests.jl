"""
This is the main test file that is called first when running any tests.
Uses modules in the test files to avoid polluting the global namespace, with QuantumTransport, Test, and other packages used. 
"""

using Test

@testset "QuantumTransport.jl Tests" begin
    println("\033[1m---------------RUNNING DRIVER TESTS---------------\033[0m")
    @testset "Driver Test" begin
        include("TestDriver.jl")
    end

    println("\033[1m----------RUNNING MATRIX INVERSION TESTS----------\033[0m")
    @testset "Matrix Inversion Test" begin
        include("TestMatrices.jl")
    end

    println("\033[1m------------RUNNING INPUT OUTPUT TESTS------------\033[0m")
    @testset "Input Output Test" begin
        include("TestIO.jl")
    end

    println("\033[1m---------RUNNING DATA VISUALIZATION TESTS----------\033[0m")
    @testset "Data Visualization Test" begin
        include("TestVisualization.jl")
    end

    println("\033[1m------------RUNNING SELF-ENERGY TESTS------------\033[0m")
    @testset "Self-Energy Test" begin
        include("TestSelfEnergies.jl")
    end

    println("\033[1m--------------RUNNING HOPPING TESTS--------------\033[0m")
    @testset "Hopping Test" begin
        include("TestHoppings.jl")
    end

end
