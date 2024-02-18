using QuantumTransport
using Test

# this file includes tests that checks for correctness for the block matrix inversion methods

include("test-matrices.jl")


# Test the RGF inverse method for correctness only
# arguments to be changed
global const args = (1000, 2, 0.2001, 1e-10, 1e-10) # Full matrix size, block size, phi, eta term, zeroThreshold term
global const testMatrix = setTestMatrix()::BlockMatrix

# Main function to drive the block matrix operations
function rgf_main()
    display(testMatrix.matrix)

    # The inverse is computed correctly for the approximatedGʳ matrix for large sizes but not fullGʳ
    # @time approximatedGʳ(0.0)
    # return

    # Call the function to test for correctness
    correct = verifyCorrectness(approximatedGʳ, fullGʳ, args[1])
    println("RGF Inverse Correctness: ", correct)

    # Call the timing function - Energy of 3.0
    timeInv(3.0)

    return correct
end

# Do the same except for only the diagonal blocks
function rgf_main_diag()
    display(testMatrix.matrix)

    # The inverse is computed correctly for the approximatedGʳ matrix for large sizes but not fullGʳ
    # @time approximatedGʳ(0.0)
    # return

    # Call the function to test for correctness
    correct = verifyCorrectness(diagApproximatedGʳ, fullGʳ, args[1])
    println("RGF Inverse Correctness: ", correct)

    # Call the timing function - Energy of 3.0
    timeInv(3.0)

    return correct
end

# Test the rgf method
@test rgf_main()

# Test the method that only returns the matrix with inverted diagonal blocks
@test rgf_main_diag()

# Test the woodbury inverse method
@test block_inv_main()