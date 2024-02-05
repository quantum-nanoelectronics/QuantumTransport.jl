using SparseArrays
using LinearAlgebra

include("./TestMatrices.jl")
include("./BlockMatrix.jl")

# Create/set test matrix
function setTestMatrix()
    # This is the initially created test matrix
    # return CreateBlockMatrix(args[1], args[2], args[3], args[5])  
    # effective mass test matrix.
    # return CreateBlockMatrix(args[1], args[2], args[3], args[5], sparse(effectiveMass(args[1] ÷ args[2], args[2])))
    # spin orbit hamiltonian test matrix - Block sizes of 2s only
    return CreateBlockMatrix(args[1], args[2], args[3], args[5], sparse(spinOrbitHamiltonian(args[1] ÷ args[2])))
end

# arguments to be changed
global const args = (1000, 2, 0.2001, 1e-10, 1e-10) # Full matrix size, block size, phi, eta term, zeroThreshold term
global const testMatrix = setTestMatrix()::BlockMatrix

# Function to compute the inverse of a block matrix using the recursive Green's function (RGF) method.
function getInvRGF(matrixObject::BlockMatrix)
    # Compute the forward and backward generators.
    forwardGen, backwardGen = computeGenerators(matrixObject)

    sparseBuild = SparseBuilder(Int[], Int[], ComplexF64[])
    # Compute the diagonal blocks of the inverse matrix.
    firstDiag, lastDiag = computeDiagonalBlocks(matrixObject, backwardGen, sparseBuild)
    computeTopBlocks(matrixObject, forwardGen, firstDiag, sparseBuild)
    computeBottomBlocks(matrixObject, backwardGen, lastDiag, sparseBuild)
    # display(sparse(sparseBuild.rowIndices, sparseBuild.columnIndices, sparseBuild.values, matrixObject.matrixSize, matrixObject.matrixSize))
    
    matrix = deepcopy(matrixObject)
    matrix.matrix = sparse(sparseBuild.rowIndices, sparseBuild.columnIndices, sparseBuild.values, matrixObject.matrixSize, matrixObject.matrixSize)
    # Return the diagonal blocks along with the top and bottom row blocks of the inverse matrix.
    # return diag, computeTopBlocks(matrixObject, forwardGen, diag[1]), computeBottomBlocks(matrixObject, backwardGen, diag[end])
    return matrix
end

function getInvJulia(matrixObject::BlockMatrix)
    matrix = deepcopy(matrixObject)
    matrix.matrix = inv(Matrix(matrixObject.matrix))
    return matrix
end

function approximatedGʳ(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + args[4]) * I - testMatrix.matrix
    return getInvRGF(matrixCopy).matrix
end

function fullGʳ(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + args[4]) * I - testMatrix.matrix
    return getInvJulia(matrixCopy).matrix
end


# Function to time the two methods of computing the inverse of a block matrix.
function timeInv(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + args[4]) * I - testMatrix.matrix
    # Time the computation of the inverse of the dense matrix using built-in inversion
    print("Julia Inverse Time: ")
    @time juliaInv = getInvJulia(matrixCopy)
    # Time the computation of the inverse of the block matrix using RGF method
    print("RGF Inverse Time: ")
    @time rgfInv = getInvRGF(matrixCopy)
end


function debugAllValues()
    a = fullGʳ(0.0)[:]
    b = approximatedGʳ(0.0)[:]
    for i in 1:size(a)[1]
        println(a[i], " ", b[i])
    end
end

function debugValues()
    # Print out comparisons between the computed blocks and the corresponding blocks in the dense inverse matrix
    # Comparing Diagonals
    println("Comparing Diagonals")
    # println(length(matrix[i:(i+numBlocks), i:(i+numBlocks)][:]))
    # println(length(diag[:]))
    for i in 1:matrixObject.numBlocks
        start_index = (i - 1) * matrixObject.blockSize + 1
        end_index = i * matrixObject.blockSize
        println(matrix[start_index:end_index, start_index:end_index])
        println(diag[i][:])
        println()
    end
    println()
    println()
    println()
    
    # Comparing Top Row
    println("Comparing Top")
    for i in 1:matrixObject.numBlocks
        start_index = (i - 1) * matrixObject.blockSize + 1
        end_index = i * matrixObject.blockSize
        println(matrix[1:matrixObject.blockSize, start_index:end_index])
        println(top[i][:])
        println()
    end
    println()
    println()
    println()
    
    # Comparing Bottom Row
    println("Comparing Bottom")
    for i in 1:matrixObject.numBlocks
        start_index = (i - 1) * matrixObject.blockSize + 1
        end_index = i * matrixObject.blockSize
        println(matrix[end-matrixObject.blockSize+1:end, start_index:end_index])
        println(bottom[i][:])
        println()
    end
    println()
    println()
    println()
end

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

