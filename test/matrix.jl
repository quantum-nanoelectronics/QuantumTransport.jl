# this file includes tests that checks for correctness for the block matrix inversion methods

# Include the necessary files 
include("test-matrices.jl")

"""
    rgf_main()

This function is the main entry point for the RGF (Recursive Green's Function) algorithm.
"""
function rgf_main(argsMatrix, testMatrix)
    display(testMatrix.matrix)

    # The inverse is computed correctly for the approximatedGʳ matrix for large sizes but not fullGʳ
    # @time approximatedGʳ(0.0)
    # return

    # Call the function to test for correctness
    correct = verifyCorrectness(approximatedGʳ, fullGʳ, argsMatrix[1])
    println("RGF Inverse Correctness: ", correct)

    # Call the timing function - Energy of 3.0
    timeInv(3.0)

    return correct
end

"""
rgf_main_diag()

This function calculates the main diagonal elements only of the RGF matrix.
"""
function rgf_main_diag(argsMatrix, testMatrix)
    display(testMatrix.matrix)

    # The inverse is computed correctly for the approximatedGʳ matrix for large sizes but not fullGʳ
    # @time approximatedGʳ(0.0)
    # return

    # Call the function to test for correctness
    correct = verifyCorrectness(diagApproximatedGʳ, fullGʳ, argsMatrix[1])
    println("RGF Inverse Correctness: ", correct)

    # Call the timing function - Energy of 3.0
    timeInv(3.0)

    return correct
end


"""
block_inv_main()

This function is the main entry point for the woodbury block inversion algorithm.
"""
function block_inv_main()
    matrix = generate_matrix()

    print("Time to compute the inverse using Julia's inv() function: ")
    @time juliaInv = inv(matrix)

    print("Time to compute the inverse using the woodbury block inversion algorithm: ")
    @time blockInv = block_inversion(matrix)

    # Calculate the norm of the difference between the two inverses
    norm_diff = norm(juliaInv - blockInv)

    # Check if the norm is close to zero when rounded to an integer
    println(round(Int, abs(norm_diff)) == 0 ? "Accurate Output." : "Inaccurate Output.")
    return round(Int, abs(norm_diff)) == 0

end


function runMatrixTests()
    # Test the RGF inverse method for correctness only
    # arguments to be changed
    # Full matrix size, block size, phi, eta term, zeroThreshold term, σ₂
    
    for matrixIndex in 0:2
        argsMatrix = (1000, 2, 0.2001, 1e-10, 1e-10, [0 -im; im 0], matrixIndex) 
        testMatrix = setVars(argsMatrix)

        println("Variables assigned")

        # Test the rgf method
        @test rgf_main(argsMatrix, testMatrix)

        # Test the method that only returns the matrix with inverted diagonal blocks
        @test rgf_main_diag(argsMatrix, testMatrix)
    end

    # Test the woodbury inverse method
    @test block_inv_main()

end

runMatrixTests()