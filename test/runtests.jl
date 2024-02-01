using QuantumTransport
using Test


# This is the main test file that is called first when running any tests.

# Include other test files here, organizing them as needed
# If you have more test files, include them in the same way:
# include("another_test_file.jl")

# Main test set that wraps all included tests
@testset "QuantumTransport.jl Tests" begin

    # This will include and run the tests from `hello_world.jl`
    println("---------------RUNNING SAMPLE TESTS---------------")
    @testset "Sample Test" begin
        include("hello_world.jl")
        include("column_major.jl")
    end
    println("----------RUNNING MATRIX INVERSION TESTS----------")
    @testset "Matrix Inversion Test" begin
        include("matrix.jl")
    end

    println("------------RUNNING INPUT OUTPUT TESTS------------")
    @testset "Input Output Test" begin
        include("io.jl")
    end

    # cannot test Data Visualization because GitHub does not have a GPU
    # @testset "Data Visualization Test" begin
    #     include("visualization.jl")
    # end

end
