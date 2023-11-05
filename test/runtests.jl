using QuantumTransport
using Test

# This is the main test file that is called first when running any tests.

# Include other test files here, organizing them as needed
# If you have more test files, include them in the same way:
# include("another_test_file.jl")

# Main test set that wraps all included tests
@testset "QuantumTransport.jl Tests" begin
    # This will include and run the tests from `hello_world.jl`
    @testset "Sample Test" begin
        include("hello_world.jl")
    end

    @testset "Column Major Test" begin
        include("TEST_column_major.jl")
    end

    # If you have other testsets, you can include them similarly:
    # @testset "Advanced Functionality" begin
    #     include("advanced_tests.jl")
    # end

end
