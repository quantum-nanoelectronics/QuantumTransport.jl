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
        include("hello-world.jl")
        include("column-major.jl")
    end
    println("----------RUNNING MATRIX INVERSION TESTS----------")
    @testset "Matrix Inversion Test" begin
        include("matrix.jl")
    end

    println("------------RUNNING INPUT OUTPUT TESTS------------")
    @testset "Input Output Test" begin
        include("io.jl")
    end

    # cannot fully test Data Visualization because GitHub does not have a GPU
    # non interactive image plots generated, however we want to generate interactive plots
    @testset "Data Visualization Test" begin
        include("visualization.jl")
    end

end
