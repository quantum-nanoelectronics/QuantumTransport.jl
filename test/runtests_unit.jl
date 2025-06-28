using Test

@testset "QuantumTransport.jl Unit Tests" begin

    println("\033[1m----------RUNNING MATRIX INVERSION TESTS----------\033[0m")
    @testset "Matrix Inversion Test" begin
        include("TestMatrices.jl")
    end

    println("\033[1m------------RUNNING INPUT OUTPUT TESTS------------\033[0m")
    @testset "Input Output Test" begin
        include("TestIO.jl")
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
