using SparseArrays
using LinearAlgebra

# Definition of a mutable structure to represent a block matrix
mutable struct SparseBlockMatrix
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

function +(a::SparseBlockMatrix, b::AbstractMatrix)
    a.matrix = a.matrix + b
    return a
end

function -(a::SparseBlockMatrix, b::AbstractMatrix)
    a.matrix = a.matrix - b
    return a
end

function -(b::AbstractMatrix, a::SparseBlockMatrix)
    a.matrix = b - a.matrix
    return a
end

# Converts to SparseBlockMatrix
function ToSparseBlockMatrix(matrix::SparseMatrixCSC{Complex{Float64},Int}, n::Int, blockSize::Int)
    if n % blockSize != 0
        error("Matrix size n must be a multiple of block size.")
    end
    numBlocks = n ÷ blockSize  # Calculate number of blocks

    return SparseBlockMatrix(matrix, n, blockSize, numBlocks)
end

# Constructor for a symmetric, block tridiagonal matrix. Returns a new SparseBlockMatrix.
function CreateSparseBlockMatrix(n::Int, blockSize::Int, phi::Float64)
    # Ensure the matrix size is a multiple of the block size
    if n % blockSize != 0
        error("Matrix size n must be a multiple of block size.")
    end
    numBlocks = n ÷ blockSize  # Calculate number of blocks

    ### The code below creates a block tridiagonal Matrix

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
    return SparseBlockMatrix(S, n, blockSize, numBlocks)
end

# Decomposes a matrix into its diagonal, top, and bottom blocks
function decomposeMatrixOld(matrixObject::SparseBlockMatrix, str::String, i::Int, matrix::Matrix{ComplexF64}, sparseBuild::SparseBuilder)
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

# TESTING - better performance than original
# Decomposes a matrix into its diagonal, top, and bottom blocks
function decomposeMatrix(matrixObject::SparseBlockMatrix, str::String, n::Int, matrix::Matrix{ComplexF64}, sparseBuild::SparseBuilder)
    threshold = 1e-10
    max_elements = length(matrix)
    rowIndices = Vector{Int}(undef, max_elements)
    colIndices = Vector{Int}(undef, max_elements)
    values = Vector{ComplexF64}(undef, max_elements)

    counter = 0
    for j in 1:size(matrix, 2)  
        for i in 1:size(matrix, 1) 
            value = matrix[j, i]  # Done intentionally for transpose, adjusted indices
            if abs(value) > threshold
                counter += 1
                rowIndices[counter] = i
                colIndices[counter] = j
                values[counter] = value
            end
        end
    end

    resize!(rowIndices, counter)
    resize!(colIndices, counter)
    resize!(values, counter)

    offset = (n - 1) * matrixObject.blockSize
    colIndices .+= offset
    if str == "diagonal"
        rowIndices .+= offset
    elseif str == "bottom"
        rowIndices .+= matrixObject.matrixSize - matrixObject.blockSize
    elseif str != "top"
        error("Invalid string")
    end

    append!(sparseBuild.rowIndices, rowIndices)
    append!(sparseBuild.columnIndices, colIndices)
    append!(sparseBuild.values, values)
end

# Retrieves the i-th diagonal block from the block matrix.
function getIthDiagonalBlock(matrixObject::SparseBlockMatrix, i::Int)
    # Calculate the start and end row indices for the i-th diagonal block
    startRow = (i - 1) * matrixObject.blockSize + 1
    endRow = i * matrixObject.blockSize
    # Extract and return the i-th diagonal block as a dense matrix
    return Matrix(matrixObject.matrix[startRow:endRow, startRow:endRow])
end

# Retrieves the block above the i-th diagonal block (i-th top diagonal block).
function getIthTopDiagonalBlock(matrixObject::SparseBlockMatrix, i::Int)
    # Do this if you want to use the bottom diagonals as the bottom ones are appromately the conjugate transpose of the top ones (testing only)
    # To always use this, this full file should be changed accordingly.
    # return getIthBottomDiagonalBlock(matrixObject, i - 1)'

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

# This function is not used. Instead, the top box is conjugate transposed and used in its place. See getIthTopDiagonalBlock().
# Retrieves the block below the i-th diagonal block (i-th bottom diagonal block).
function getIthBottomDiagonalBlock(matrixObject::SparseBlockMatrix, i::Int)
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
function getIthTopRowBlock(matrixObject::SparseBlockMatrix, i::Int)
    # Calculate the column indices for the top row blocks at the i-th position
    startCol = (i - 1) * matrixObject.blockSize + 1
    endCol = i * matrixObject.blockSize
    # The top row block count is equal to the block size
    topRowsCount = matrixObject.blockSize
    # Extract and return the top row block as a dense matrix
    return Matrix(matrixObject.matrix[1:topRowsCount, startCol:endCol])
end

# Retrieves the bottom row block at the i-th position.
function getIthBottomRowBlock(matrixObject::SparseBlockMatrix, i::Int)
    # Calculate the column indices for the bottom row blocks at the i-th position
    startCol = (i - 1) * matrixObject.blockSize + 1
    endCol = i * matrixObject.blockSize
    # Calculate the starting row index for the bottom row blocks
    bottomRowsStart = size(matrixObject.matrix, 1) - matrixObject.blockSize + 1
    # Extract and return the bottom row block as a dense matrix
    return Matrix(matrixObject.matrix[bottomRowsStart:end, startCol:endCol])
end

# Function to compute forward and backward generators for a block matrix.
function computeGenerators(matrixObject::SparseBlockMatrix)
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
function computeDiagonalBlocks(matrixObject::SparseBlockMatrix, backwardGen::Vector{Matrix{ComplexF64}}, sparseBuild::SparseBuilder)
    # Initialize diagonal blocks as an array of zero matrices.
    # diagonalBlocks = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 1:matrixObject.numBlocks]
    # Compute the diagonal block separately.
    # diagonalBlocks[matrixObject.numBlocks] = backwardGen[matrixObject.numBlocks]

    lastDiagBlock = backwardGen[matrixObject.numBlocks]
    currentDiagBlock = nothing
    lastBlock = lastDiagBlock
    decomposeMatrix(matrixObject, "diagonal", matrixObject.numBlocks, lastDiagBlock, sparseBuild)
    
    for i in (matrixObject.numBlocks - 1:-1:1)
        topBlock = getIthTopDiagonalBlock(matrixObject, i + 1)
        currentDiagBlock = backwardGen[i] * (I + topBlock * lastDiagBlock * topBlock' * backwardGen[i])
        decomposeMatrix(matrixObject, "diagonal", i, currentDiagBlock, sparseBuild) 
        lastDiagBlock = currentDiagBlock
    end

    firstBlock = currentDiagBlock
    return firstBlock, lastBlock
end


# Function to compute the top row of blocks in the inverse of a block matrix.
function computeTopBlocks(matrixObject::SparseBlockMatrix, forwardGen::Vector{Matrix{ComplexF64}}, firstDiagBlock::Matrix{ComplexF64}, sparseBuild::SparseBuilder)
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
function computeBottomBlocks(matrixObject::SparseBlockMatrix, backwardGen::Vector{Matrix{ComplexF64}}, lastDiagBlock::Matrix{ComplexF64}, sparseBuild::SparseBuilder)
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
