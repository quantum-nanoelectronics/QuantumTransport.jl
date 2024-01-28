using Arpack
using SparseArrays
using LinearAlgebra
import Base.summarysize

# Definition of a mutable structure to represent a block matrix
mutable struct BlockMatrix
    matrix::SparseMatrixCSC{Complex{Float64},Int}  # The sparse block tridiagonal matrix
    matrixSize::Int  # Total size of the matrix
    blockSize::Int  # Size of each square block within the matrix
    numBlocks::Int  # Number of blocks along one dimension
end

mutable struct SparseBuilder
    rowIndices::Vector{Int}
    columnIndices::Vector{Int}
    values::Vector{ComplexF64}
end

# Constructor for a symmetric, block tridiagonal matrix
function BlockMatrix(n::Int, blockSize::Int, phi::Float64)
    # Ensure the matrix size is a multiple of the block size
    if n % blockSize != 0
        error("Matrix size n must be a multiple of block size.")
    end

    numBlocks = n ÷ blockSize  # Calculate number of blocks

    # Preallocate vectors for SparseMatrixCSC construction
    vals = Complex{Float64}[]
    rows = Int[]
    cols = Int[]

    # Loop to fill in the blocks of the matrix
    for blockIndex in 1:numBlocks
        # Determine the range of indices for the current block
        startIdx = (blockIndex - 1) * blockSize + 1
        endIdx = blockIndex * blockSize

        # Fill the diagonal block with a constant value
        for i in startIdx:endIdx
            push!(vals, 2.0 + i)  # Diagonal value
            push!(rows, i)     # Row index
            push!(cols, i)     # Column index
        end

        # Fill the off-diagonal blocks if not at the boundary
        if blockIndex < numBlocks
            for i in startIdx:endIdx
                j = startIdx + blockSize + (endIdx - i)
                push!(vals, exp(-1im * phi))  # Upper diagonal
                # push!(vals, 1im)  # Positive imaginary unit above the diagonal
                # push!(vals, 1)   # Off-diagonal value above the main diagonal
                push!(rows, i)   # Row index
                push!(cols, j)   # Column index

                push!(vals, exp(1im * phi)) # Lower diagonal
                # push!(vals, -1im)  # Negative imaginary unit below the diagonal
                # push!(vals, 1)   # Off-diagonal value below the main diagonal
                push!(rows, j)   # Column index
                push!(cols, i)   # Row index
            end
        end
    end

    # Construct the sparse matrix from the populated vectors
    S = sparse(rows, cols, vals, n, n)
    return BlockMatrix(S, n, blockSize, numBlocks)
end

function decomposeMatrix(matrixObject::BlockMatrix, str::String, i::Int, matrix::Matrix{ComplexF64}, sparseBuild::SparseBuilder)
    threshold = 1e-10  # Define a small threshold
    rowIndices = Vector{Int}()
    colIndices = Vector{Int}()
    values = Vector{ComplexF64}()

    for col in 1:size(matrix, 2)
        for row in 1:size(matrix, 1)
            value = matrix[col, row] # this is done on purpose for the transpose
            if abs(value) > threshold  # Check if the magnitude of the value is above the threshold
                push!(rowIndices, row)
                push!(colIndices, col)
                push!(values, value)
            end
        end
    end

    offset = (i - 1) * matrixObject.blockSize
    colIndices .+= offset
    if str == "diagonal"
        rowIndices .+= offset
    elseif str == "top"
        # No change to rowIndices
        nothing
    elseif str == "bottom"
        rowIndices .+= matrixObject.matrixSize - matrixObject.blockSize
    else
        error("Invalid string")
    end

    append!(sparseBuild.rowIndices, rowIndices)
    append!(sparseBuild.columnIndices, colIndices)
    append!(sparseBuild.values, values)
end



# This function is not used. Instead, the top box is conjugate transposed and used in its place.
# Retrieves the i-th diagonal block from the block matrix.
function getIthDiagonalBlock(matrixObject::BlockMatrix, i::Int)
    # Calculate the start and end row indices for the i-th diagonal block
    startRow = (i - 1) * matrixObject.blockSize + 1
    endRow = i * matrixObject.blockSize
    # Extract and return the i-th diagonal block as a dense matrix
    return Matrix(matrixObject.matrix[startRow:endRow, startRow:endRow])
end

# Retrieves the block above the i-th diagonal block (i-th top diagonal block).
function getIthTopDiagonalBlock(matrixObject::BlockMatrix, i::Int)
    if i == 1
        # The first block does not have a top diagonal block
        error("There is no top block for the first block")
    end
    # Calculate the row and column indices for the i-th top diagonal block
    startRow = (i - 1) * matrixObject.blockSize + 1
    endRow = i * matrixObject.blockSize
    startCol = (i - 2) * matrixObject.blockSize + 1
    endCol = (i - 1) * matrixObject.blockSize
    # Extract and return the i-th top diagonal block as a dense matrix
    return Matrix(matrixObject.matrix[startRow:endRow, startCol:endCol])
end

# Retrieves the block below the i-th diagonal block (i-th bottom diagonal block).
function getIthBottomDiagonalBlock(matrixObject::BlockMatrix, i::Int)
    if i == matrixObject.numBlocks
        # The last block does not have a bottom diagonal block
        error("There is no bottom block for the last block")
    end
    # Calculate the row and column indices for the i-th bottom diagonal block
    startRow = i * matrixObject.blockSize + 1
    endRow = (i + 1) * matrixObject.blockSize
    startCol = (i - 1) * matrixObject.blockSize + 1
    endCol = i * matrixObject.blockSize
    # Extract and return the i-th bottom diagonal block as a dense matrix
    return Matrix(matrixObject.matrix[startRow:endRow, startCol:endCol])
end

# Retrieves the top row block at the i-th position.
function getIthTopRowBlock(matrixObject::BlockMatrix, i::Int)
    # Calculate the column indices for the top row blocks at the i-th position
    startCol = (i - 1) * matrixObject.blockSize + 1
    endCol = i * matrixObject.blockSize
    # The top row block count is equal to the block size
    topRowsCount = matrixObject.blockSize
    # Extract and return the top row block as a dense matrix
    return Matrix(matrixObject.matrix[1:topRowsCount, startCol:endCol])
end

# Retrieves the bottom row block at the i-th position.
function getIthBottomRowBlock(matrixObject::BlockMatrix, i::Int)
    # Calculate the column indices for the bottom row blocks at the i-th position
    startCol = (i - 1) * matrixObject.blockSize + 1
    endCol = i * matrixObject.blockSize
    # Calculate the starting row index for the bottom row blocks
    bottomRowsStart = size(matrixObject.matrix, 1) - matrixObject.blockSize + 1
    # Extract and return the bottom row block as a dense matrix
    return Matrix(matrixObject.matrix[bottomRowsStart:end, startCol:endCol])
end


# Function to compute forward and backward generators for a block matrix.
function computeGenerators(matrixObject::BlockMatrix)
    # Initialize forward and backward generators as arrays of zero matrices.
    forwardGen = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 1:matrixObject.numBlocks]
    backwardGen = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 1:matrixObject.numBlocks]
    # Initialize the first backward generator using the first block's inverse and the top diagonal block.
    backwardGen[1] = inv(getIthDiagonalBlock(matrixObject, 1))
    # Calculate remaining backward generators using a loop.
    for i in 2:(matrixObject.numBlocks)
        backwardGen[i] = inv(getIthDiagonalBlock(matrixObject, i) - (getIthTopDiagonalBlock(matrixObject, i) * backwardGen[i-1] * getIthTopDiagonalBlock(matrixObject, i)'))
    end 
    # Initialize the last forward generator using the last block's top diagonal block and its inverse.
    forwardGen[matrixObject.numBlocks] = inv(getIthDiagonalBlock(matrixObject, matrixObject.numBlocks))
    # Calculate remaining forward generators in reverse order using a loop.
    for i in (matrixObject.numBlocks - 1):-1:1
        forwardGen[i] = inv(getIthDiagonalBlock(matrixObject, i) - (getIthTopDiagonalBlock(matrixObject, i + 1) * forwardGen[i + 1] * getIthTopDiagonalBlock(matrixObject, i + 1)'))
    end
    return forwardGen, backwardGen
end

# Function to compute diagonal blocks of a block matrix using forward generators.
function computeDiagonalBlocks(matrixObject::BlockMatrix, backwardGen::Vector{Matrix{ComplexF64}}, sparseBuild::SparseBuilder)
    # Initialize diagonal blocks as an array of zero matrices.
    # diagonalBlocks = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 1:matrixObject.numBlocks]
    # Compute the diagonal block separately.
    # diagonalBlocks[matrixObject.numBlocks] = backwardGen[matrixObject.numBlocks]

    lastDiagBlock = backwardGen[matrixObject.numBlocks]
    currentDiagBlock = nothing
    lastBlock = lastDiagBlock
    decomposeMatrix(matrixObject, "diagonal", matrixObject.numBlocks, lastDiagBlock, sparseBuild)
    
    for i in (matrixObject.numBlocks - 1:-1:1)
        currentDiagBlock = backwardGen[i] * (I + getIthTopDiagonalBlock(matrixObject, i + 1) * lastDiagBlock * getIthTopDiagonalBlock(matrixObject, i + 1)' * backwardGen[i])
        decomposeMatrix(matrixObject, "diagonal", i, currentDiagBlock, sparseBuild) 
        lastDiagBlock = currentDiagBlock
    end

    firstBlock = currentDiagBlock
    return firstBlock, lastBlock
end

# Function to compute the top row of blocks in the inverse of a block matrix.
function computeTopBlocks(matrixObject::BlockMatrix, forwardGen::Vector{Matrix{ComplexF64}}, firstDiagBlock::Matrix{ComplexF64}, sparseBuild::SparseBuilder)
    # Initialize top row blocks as an array of zero matrices.
    # topBlocks = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 1:matrixObject.numBlocks]
    # Set the first block in the top row.

    lastTopBlock = firstDiagBlock
    for i in 2:matrixObject.numBlocks
        currentTopBlock = (forwardGen[i] * -getIthTopDiagonalBlock(matrixObject, i)' * lastTopBlock)
        # println(currentTopBlock)

        decomposeMatrix(matrixObject, "top", i, currentTopBlock, sparseBuild) 
        lastTopBlock = currentTopBlock
    end
end

# Function to compute the bottom row of blocks in the inverse of a block matrix.
function computeBottomBlocks(matrixObject::BlockMatrix, backwardGen::Vector{Matrix{ComplexF64}}, lastDiagBlock::Matrix{ComplexF64}, sparseBuild::SparseBuilder)
    # Initialize bottom row blocks as an array of zero matrices.
    # bottomBlocks = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 1:matrixObject.numBlocks]
    # Set the last block in the bottom row.
    # bottomBlocks[matrixObject.numBlocks] = lastDiagBlock
    lastBottomBlock = lastDiagBlock
    # Compute each bottom block iteratively in reverse order.
    for i in (matrixObject.numBlocks - 1):-1:1
        currentBottomBlock = backwardGen[i] * -getIthTopDiagonalBlock(matrixObject, i + 1) * lastBottomBlock
        decomposeMatrix(matrixObject, "bottom", i, currentBottomBlock, sparseBuild) 
        lastBottomBlock = currentBottomBlock
    end
    # Return the array of bottom row blocks.
    # pop!(bottomBlocks)
end

# Function to convert the block matrix to a dense matrix format. Recommended only for small matrices due to high space complexity cost.
function getDense(matrixObject::BlockMatrix)
    # Convert and return the block matrix as a dense matrix.
    return Matrix(matrixObject.matrix)
end

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

    
    # Return the diagonal blocks along with the top and bottom row blocks of the inverse matrix.
    # return diag, computeTopBlocks(matrixObject, forwardGen, diag[1]), computeBottomBlocks(matrixObject, backwardGen, diag[end])
    return BlockMatrix(sparse(sparseBuild.rowIndices, sparseBuild.columnIndices, sparseBuild.values, matrixObject.matrixSize, matrixObject.matrixSize), matrixObject.matrixSize, matrixObject.blockSize, matrixObject.numBlocks)
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
function main()
    args = (8, 2, 0.2001)
    matrixObject = BlockMatrix(args...)  # Create a block matrix of size 1000 with block size 1 0.2001
    dense = getDense(matrixObject)  # Convert the block matrix to a dense matrix

    # Time the computation of the inverse of the dense matrix using built-in inversion
    juliaInv = inv(dense)

    # Time the computation of the inverse of the block matrix using RGF method
    rgfInv = getInvRGF(matrixObject)

    return

end

main()