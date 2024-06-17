using SparseArrays
using LinearAlgebra

include("SparseBlockMatrix.jl")

function getInvJulia!(matrixObject::SparseBlockMatrix)
    matrixObject.matrix = inv(Matrix(matrixObject.matrix))
end

# Function to compute the inverse of a block matrix using the recursive Green's function (RGF) method.
function getInvRGF!(matrixObject::SparseBlockMatrix)
    # does not work if there is only one block
    if matrixObject.numBlocks == 1
        getInvJulia!(matrixObject)
        return
    end

    # Compute the forward and backward generators.
    forwardGen, backwardGen = computeGenerators(matrixObject)

    sparseBuild = SparseBuilder(Int[], Int[], ComplexF64[])
    # Compute the diagonal blocks of the inverse matrix.
    firstDiag, lastDiag = computeDiagonalBlocks(matrixObject, backwardGen, sparseBuild)
    computeTopBlocks(matrixObject, forwardGen, firstDiag, sparseBuild)
    computeBottomBlocks(matrixObject, backwardGen, lastDiag, sparseBuild)
    # display(sparse(sparseBuild.rowIndices, sparseBuild.columnIndices, sparseBuild.values, matrixObject.matrixSize, matrixObject.matrixSize))
    
    matrixObject.matrix = sparse(sparseBuild.rowIndices, sparseBuild.columnIndices, sparseBuild.values, matrixObject.matrixSize, matrixObject.matrixSize)
    # Return the diagonal blocks along with the top and bottom row blocks of the inverse matrix.
    # return diag, computeTopBlocks(matrixObject, forwardGen, diag[1]), computeBottomBlocks(matrixObject, backwardGen, diag[end])
end

# Function to perform the getInvRGF function except without the top and bottom rows of the matrix
function getInvRGFDiagonal!(matrixObject::SparseBlockMatrix)
    # does not work if there is only one block
    if matrixObject.numBlocks == 1
        getInvJulia!(matrixObject)
        return
    end

    # Compute the forward and backward generators.
    _, backwardGen = computeGenerators(matrixObject)

    sparseBuild = SparseBuilder(Int[], Int[], ComplexF64[])
    # Compute the diagonal blocks of the inverse matrix.
    computeDiagonalBlocks(matrixObject, backwardGen, sparseBuild)
    
    matrixObject.matrix = sparse(sparseBuild.rowIndices, sparseBuild.columnIndices, sparseBuild.values, matrixObject.matrixSize, matrixObject.matrixSize)
end


function getInvRGF(matrixObject::SparseBlockMatrix)
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

# Function to perform the getInvRGF function except without the top and bottom rows of the matrix
function getInvRGFDiagonal(matrixObject::SparseBlockMatrix)
    # Compute the forward and backward generators.
    _, backwardGen = computeGenerators(matrixObject)

    sparseBuild = SparseBuilder(Int[], Int[], ComplexF64[])
    # Compute the diagonal blocks of the inverse matrix.
    computeDiagonalBlocks(matrixObject, backwardGen, sparseBuild)
    
    matrix = deepcopy(matrixObject)
    matrix.matrix = sparse(sparseBuild.rowIndices, sparseBuild.columnIndices, sparseBuild.values, matrixObject.matrixSize, matrixObject.matrixSize)
    return matrix
end

function getInvJulia(matrixObject::SparseBlockMatrix)
    matrix = deepcopy(matrixObject)
    matrix.matrix = inv(Matrix(matrixObject.matrix))
    return matrix
end