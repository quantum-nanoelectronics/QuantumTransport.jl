# this file includes tests that checks for correctness for the block matrix inversion methods

module TestMatrices
# Usual imports
using QuantumTransport
using Test

# Other imports
using LinearAlgebra
using SparseArrays
using BenchmarkTools

# This data is to be set in this code
# Included in the module which will not pollute the global namespace
argsMatrix = nothing
testMatrix = nothing

# Include the necessary files 
include("TestMatricesHelper.jl")


"""
rgf_main(argsMatrix, testMatrix, diag=false)

The main function for testing the RGF algorithm for correctness and timing.

## Arguments
- `argsMatrix`: A tuple containing the arguments for the RGF matrix algorithm.
- `testMatrix`: A matrix used for testing purposes.
- `diag`: A boolean indicating whether to compute the diagonal elements of the inverse matrix. Default is `false`.

## Returns
- `correct`: A boolean indicating whether the RGF inverse is correct.

"""
function rgf_main(argsMatrix, testMatrix, diag::Bool=false)
    println("\033[1mTesting the RGF algorithm for correctness and timing. Diagonal only? : $diag \033[0m")
    display(testMatrix.matrix)

    # The inverse is computed correctly for the approximatedGʳ matrix for large sizes but not fullGʳ

    # Call the function to test for correctness
    if diag
        correct = verifyCorrectness(diagApproximatedGʳ, fullGʳ, argsMatrix[1])
    else
        correct = verifyCorrectness(approximatedGʳ, fullGʳ, argsMatrix[1])
    end
    println("RGF Inverse Correctness: ", correct)

    # Call the timing function at an energy level - argsMatrix[7]
    @test timeInv(argsMatrix[7])

    println()

    return correct
end


"""
block_inv_main()

This function is the main entry point for the woodbury block inversion algorithm.
"""
function block_inv_main()
    println("\033[1mTesting the woodbury block inversion algorithm for correctness\033[0m")

    matrix = generate_matrix()

    print("Time to compute the inverse using Julia's inv() function: ")
    @time juliaInv = inv(matrix)

    print("Time to compute the inverse using the woodbury block inversion algorithm: ")
    @time blockInv = block_inversion(matrix)

    # Calculate the norm of the difference between the two inverses
    norm_diff = norm(juliaInv - blockInv)

    # Check if the norm is close to zero when rounded to an integer
    println((Int, abs(norm_diff))) # == 0 ? "Accurate Output." : "Inaccurate Output.")
    return round(Int, abs(norm_diff)) == 0

end


"""
    runMatrixTests()

Run a series of tests on matrix operations.

This function tests the RGF (Recursive Green's Function) inverse method for correctness and the woodbury inverse method for correctness.

# Arguments
- `matrixIndex::Int`: The index of the matrix to be tested.

"""
function runMatrixTests()
    # Test the RGF inverse method for correctness only
    # arguments to be changed
    # Full matrix size, block size, phi, eta term, zeroThreshold term, σ₂, energy, 
    # ϕ
    # η

    ### 
    # Commented for now for time purposes
    # for matrixIndex in 0:2
    #     global argsMatrix = (1000, 2, 0.2001, 1e-10, 1e-10, [0 -im; im 0], 3.0, matrixIndex)
    #     global testMatrix = setVars(argsMatrix)

    #     # Test the rgf method
    #     @test rgf_main(argsMatrix, testMatrix)

    #     # Test the method that only returns the matrix with inverted diagonal blocks
    #     # @test rgf_main(argsMatrix, testMatrix, true)
    # end
    # # Test the woodbury inverse method
    # @test block_inv_main()
    ###

    global argsMatrix = (200, 2, 0.2001, 1e-10, 1e-10, [0 -im; im 0], 3.0, 1)
    global testMatrix = setVars(argsMatrix)
    @test rgf_main(argsMatrix, testMatrix)

end

runMatrixTests()

end # module TestMatrices
