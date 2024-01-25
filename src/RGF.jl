using Arpack
using SparseArrays
using LinearAlgebra
import Base.summarysize

mutable struct BlockMatrix
    matrix::SparseMatrixCSC{Complex{Float64},Int}
    matrixSize::Int
    blockSize::Int
    numBlocks::Int
end

# symmetric, block tri-diagonal matrix
function BlockMatrix(n::Int, blockSize::Int)
    # Ensure the matrix size is a multiple of the block size
    if n % blockSize != 0
        error("Matrix size n must be a multiple of block size.")
    end

    # Number of blocks along one dimension
    numBlocks = n ÷ blockSize

    # Preallocate vectors for SparseMatrixCSC
    vals = Complex{Float64}[]
    rows = Int[]
    cols = Int[]

    for blockIndex in 1:numBlocks
        # The starting and ending index for the current block
        startIdx = (blockIndex - 1) * blockSize + 1
        endIdx = blockIndex * blockSize

        # Fill the diagonal block
        for i in startIdx:endIdx
            push!(vals, -2.0)
            push!(rows, i)
            push!(cols, i)
        end

        # Fill the off-diagonal blocks if not at the boundary
        if blockIndex < numBlocks
            for i in startIdx:endIdx
                j = startIdx + blockSize + (endIdx - i)
                # push!(vals, 1im)  # Positive imaginary unit above the diagonal
                push!(vals, 1)
                push!(rows, i)
                push!(cols, j)

                # push!(vals, -1im)  # Negative imaginary unit below the diagonal
                push!(vals, 1)
                push!(rows, j)
                push!(cols, i)
            end
        end
    end

    # Create the sparse matrix
    S = sparse(rows, cols, vals, n, n)
    return BlockMatrix(S, n, blockSize, numBlocks)
end

function getIthDiagonalBlock(matrixObject::BlockMatrix, i::Int)
    start_row = (i - 1) * matrixObject.blockSize + 1
    end_row = i *  matrixObject.blockSize
    return Matrix(matrixObject.matrix[start_row:end_row, start_row:end_row])
end

function getIthTopDiagonalBlock(matrixObject::BlockMatrix, i::Int)
    # there is no top block for the first diagonal block
    if i == 1
        error("There is no top block for the first block")
    end

    start_row = (i - 1) * matrixObject.blockSize + 1
    end_row = i * matrixObject.blockSize
    start_col = (i - 2) * matrixObject.blockSize + 1
    end_col = (i - 1) * matrixObject.blockSize

    return Matrix(matrixObject.matrix[start_row:end_row, start_col:end_col])
end

function getIthBottomDiagonalBlock(matrixObject::BlockMatrix, i::Int)
    # there is no bottom block for the last diagonal block
    if i == matrixObject.numBlocks
        error("There is no bottom block for the last block")
    end

    start_row = i * matrixObject.blockSize + 1
    end_row = (i + 1) * matrixObject.blockSize
    start_col = (i - 1) * matrixObject.blockSize + 1
    end_col = i * matrixObject.blockSize

    return Matrix(matrixObject.matrix[start_row:end_row, start_col:end_col])
end

function getIthTopRowBlock(matrixObject::BlockMatrix, i::Int)
    start_col = (i - 1) * matrixObject.blockSize + 1
    end_col = i * matrixObject.blockSize
    top_rows_count = matrixObject.blockSize
    return Matrix(matrixObject.matrix[1:top_rows_count, start_col:end_col])
end

function getIthBottomRowBlock(matrixObject::BlockMatrix, i::Int)
    start_col = (i - 1) * matrixObject.blockSize + 1
    end_col = i * matrixObject.blockSize
    bottom_rows_start = size(matrixObject.matrix, 1) - matrixObject.blockSize + 1
    return Matrix(matrixObject.matrix[bottom_rows_start:end, start_col:end_col])
end

function computeGenerators(matrixObject::BlockMatrix)
    forwardGen = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 2:matrixObject.numBlocks]
    backwardGen = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 2:matrixObject.numBlocks]

    # Filling in the backward generator matrix list
    backwardGen[1] = inv(getIthDiagonalBlock(matrixObject, 1)) * -getIthTopDiagonalBlock(matrixObject, 2)
    for i in 2:(matrixObject.numBlocks - 1)
        backwardGen[i] = inv(getIthDiagonalBlock(matrixObject, i) - (-getIthBottomDiagonalBlock(matrixObject, i - 1) * backwardGen[i-1])) * -getIthTopDiagonalBlock(matrixObject, i + 1)
    end 

    # Filling in the forward generator matrix list
    forwardGen[matrixObject.numBlocks - 1] = -getIthTopDiagonalBlock(matrixObject, matrixObject.numBlocks) * inv(getIthDiagonalBlock(matrixObject, matrixObject.numBlocks))
    for i in (matrixObject.numBlocks - 2):-1:1
        forwardGen[i] = -getIthTopDiagonalBlock(matrixObject, i + 1) * inv(getIthDiagonalBlock(matrixObject, i + 1) - (forwardGen[i + 1] * -getIthBottomDiagonalBlock(matrixObject, i + 1)))
    end

    return forwardGen, backwardGen
end


function computeDiagonalBlocks(matrixObject::BlockMatrix, forwardGen::Vector{Matrix{ComplexF64}})
    diagonalBlocks = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 1:matrixObject.numBlocks]
    diagonalBlocks[1] = inv(getIthDiagonalBlock(matrixObject, 1) - (forwardGen[1] * -getIthBottomDiagonalBlock(matrixObject, 1)))
    for i in (2:matrixObject.numBlocks - 1)
        term1 = inv(getIthDiagonalBlock(matrixObject, i) - (forwardGen[i] * -getIthBottomDiagonalBlock(matrixObject, i)))
        term2 = I + -getIthBottomDiagonalBlock(matrixObject, i - 1) * diagonalBlocks[i - 1] * forwardGen[i - 1]
        diagonalBlocks[i] = term1 * term2
    end
    diagonalBlocks[matrixObject.numBlocks] = inv(getIthDiagonalBlock(matrixObject, matrixObject.numBlocks)) * (I + (-getIthBottomDiagonalBlock(matrixObject, matrixObject.numBlocks - 1) * diagonalBlocks[matrixObject.numBlocks - 1] * forwardGen[end]))
    return diagonalBlocks
end

# Compute the top row inverse
function computeTopBlocks(matrixObject::BlockMatrix, forwardGen::Vector{Matrix{ComplexF64}}, firstDiagBlock::Matrix{ComplexF64})
    topBlocks = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 1:matrixObject.numBlocks]
    topBlocks[1] = firstDiagBlock
    runningProduct = I
    for i in 2:matrixObject.numBlocks
        runningProduct = runningProduct * forwardGen[i - 1]
        topBlocks[i] = firstDiagBlock * runningProduct
    end
    return topBlocks
end

# Compute the bottom row inverse
function computeBottomBlocks(matrixObject::BlockMatrix, backwardGen::Vector{Matrix{ComplexF64}}, lastDiagBlock::Matrix{ComplexF64})
    bottomBlocks = [zeros(Complex{Float64}, matrixObject.blockSize, matrixObject.blockSize) for i in 1:matrixObject.numBlocks]
    bottomBlocks[matrixObject.numBlocks] = lastDiagBlock
    runningProduct = I

    for i in (matrixObject.numBlocks - 1):-1:1
        runningProduct = runningProduct * backwardGen[i]
        bottomBlocks[i] = runningProduct * lastDiagBlock
    end

    # println(backwardGen[matrixObject.numBlocks - 1] * allDiagBlocks[matrixObject.numBlocks - 1])

    return bottomBlocks

end

# should only be used for smaller test cases
function getDense(matrixObject::BlockMatrix)
    return Matrix(matrixObject.matrix)
end

function getInvRGF(matrixObject::BlockMatrix)
    forwardGen, backwardGen = computeGenerators(matrixObject)
    diag = computeDiagonalBlocks(matrixObject, forwardGen)
    return diag, computeTopBlocks(matrixObject, forwardGen, diag[1]), computeBottomBlocks(matrixObject, backwardGen, diag[end])
end

function main()
    matrixObject = BlockMatrix(1000, 1)
    dense = getDense(matrixObject)

    @time diag, top, bottom = getInvRGF(matrixObject)
    
    @time matrix = inv(dense)
    rows, _ = size(matrix)


    # forwardGen, backwardGen = computeGenerators(matrixObject)


    ### RANDOM TESTS
    # println(getIthTopRowBlock(matrixObject,1))
    # println(getIthTopDiagonalBlock(matrixObject,6))
    # println(matrixObject.numBlocks)
    # println(dense)
    # println(forwardGen)
    # backwardGen[1][1] = 0.5
    # backwardGen[2][1] = 0.4
    # backwardGen[3][1] = 0.3
    # println(backwardGen)

    # return

    ### COMPARING DIAGONALS
    println("Comparing Diagonals")
    for i in 1:rows
        print(matrix[i, i])
        print("   ")
        print(diag[i])
        println()
    end
    println()
    println()
    println()

    ### COMPARING TOP ROW
    println("Comparing Top")
    for i in 1:rows
        print(matrix[1, i])
        print("   ")
        print(top[i])
        println()
    end
    println()
    println()
    println()
    
    ### COMPARING BOTTOM ROW
    println("Comparing Bottom")
    for i in 1:rows
        print(matrix[rows, i])
        print("   ")
        print(bottom[i])
        println()
    end
    println()
    println()
    println()
end


main()





































###OLD CODE BELOW














    
#calculate diagonal, top, and bottom


# function main()
#     println("Efficient RGF Method running")

#     n = 100  # Total size of the matrix
#     block_size = 5  # Size of each block
#     num_blocks = n ÷ block_size  # Number of blocks

#     # Create a large sparse tridiagonal matrix
#     sparseMatrix = spdiagm(0 => fill(1.0 + 0im, n), 1 => fill(-2.0 + 0im, n-1), -1 => fill(-2.0 + 0im, n-1))

#     # Apply efficient RGF method
#     top_blocks, bottom_blocks, diagonal_blocks = rgf_efficient(sparseMatrix, block_size, num_blocks)

#     # You can add more code here to analyze or use the parts of the Green's Function
#     # ...

# end



####
# OLD code to test size, type, shape, etc.
####

   # value = sparseMatrix[1, 1]
    # println(value)

    # i, j, k, l = 10, 15, 10, 15  
    # subBlockSparse = sparseMatrix[i:j, k:l]
    # subBlockDense = Array(subBlockSparse)
    # print(size(sparseMatrix, 1))

    # greenFunction = rgf(sparseMatrix)


function testing()
    println("RGF Method running")

    n = 1000

    vals = Complex{Float64}[]
    rows = Int[]
    cols = Int[]

    for i in 1:n
        # Main diagonal
        push!(vals, 1.0)
        push!(rows, i)
        push!(cols, i)

        # top diagonal
        if i > 1
            push!(vals, -2.0)
            push!(rows, i)
            push!(cols, i-1)
        end

        # bottom diagonal
        if i < n
            push!(vals, -2.0)
            push!(rows, i)
            push!(cols, i+1)
        end
    end

    sparseMatrix = sparse(rows, cols, vals, n, n)
    value = sparseMatrix[1, 2]



    testMatrix = zeros(Complex{Float64}, n, n)
    for i in 1:n
        testMatrix[i, i] = 1.0
        if i > 1
            testMatrix[i, i-1] = -2.0
        end
        if i < n
            testMatrix[i, i+1] = -2.0
        end
    end

    # testMatrix = [1 -2 0 0 0; -2 1 -2 0 0; 0 -2 1 -2 0; 0 0 -2 1 -2; 0 0 0 -2 1]

    complexTestMatrix = Complex{Float64}.(testMatrix)
    sparseMatrix = sparse(complexTestMatrix)
    denseConvertedMatrix = Matrix(sparseMatrix)

    println("Type of testMatrix: ", typeof(testMatrix))
    println("Type of complexTestMatrix: ", typeof(complexTestMatrix))
    println("Type of sparseMatrix: ", typeof(sparseMatrix))
    println("Type of denseConvertedMatrix: ", typeof(denseConvertedMatrix))
    println()

    println("Size of sparseMatrix (bytes): ", summarysize(sparseMatrix))
    println("Size of denseConvertedMatrix (bytes): ", summarysize(denseConvertedMatrix))

end



# -2  0
# 0 -2

# top 
# 0 -i
# i 0

# bottom 
# 0 I
# -i 0


function BlockMatrix2(n::Int, blockSize::Int)
    # Ensure the matrix size is a multiple of the block size
    if n % blockSize != 0
        error("Matrix size n must be a multiple of block size.")
    end

    # Number of blocks along one dimension
    numBlocks = n ÷ blockSize

    # Preallocate vectors for SparseMatrixCSC
    vals = Complex{Float64}[]
    rows = Int[]
    cols = Int[]

    for blockIndex in 1:numBlocks
        # The starting and ending index for the current block
        startIdx = (blockIndex - 1) * blockSize + 1
        endIdx = blockIndex * blockSize

        # Fill the diagonal block
        for i in startIdx:endIdx, j in startIdx:endIdx
            push!(vals, 1.0)
            push!(rows, i)
            push!(cols, j)
        end

        # Fill the off-diagonal blocks if not at the boundary
        if blockIndex < numBlocks
            for i in startIdx:endIdx, j in (endIdx + 1):(endIdx + blockSize)
                push!(vals, -2.0)
                push!(rows, i)
                push!(cols, j)
                push!(vals, -2.0)
                push!(rows, j)
                push!(cols, i)
            end
        end
    end

    # Create the sparse matrix
    S = sparse(rows, cols, vals, n, n)
    return BlockMatrix(S, n, blockSize, n ÷ blockSize)
end

function BlockMatrix1(n::Int, blockSize::Int)
    println("Creating test matrix")
    vals = Complex{Float64}[]
    rows = Int[]
    cols = Int[]
    for i in 1:n
        push!(vals, 1.0)
        push!(rows, i)
        push!(cols, i)
        if i > 1
            push!(vals, -2.0)
            push!(rows, i)
            push!(cols, i-1)
        end
        if i < n
            push!(vals, -2.0)
            push!(rows, i)
            push!(cols, i+1)
        end
    end
    return BlockMatrix(sparse(rows, cols, vals, n, n), n, blockSize, n ÷ blockSize)
end