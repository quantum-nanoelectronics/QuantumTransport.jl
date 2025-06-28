"""
This script runs integration tests for the QuantumTransport.jl package.
"""

using Test

@testset "QuantumTransport.jl Integration Tests" begin
    
    println("\033[1m---------------RUNNING DRIVER---------------\033[0m")
    @testset "Driver Test" begin
        include("TestDriver.jl")
    end

    println("\033[1m---------RUNNING DATA VISUALIZATION----------\033[0m")
    @testset "Data Visualization Test" begin
        include("TestVisualization.jl")
    end

end
