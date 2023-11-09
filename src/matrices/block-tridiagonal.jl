# Define a BlockTridiagonal type to hold the three diagonals
# generic type T can hold any block BlockTridiagonal matrix
# contains one diagonal sparse matrix of size n, 2 sub/super of size n-1

mutable struct BlockTridiagonal{T} <: AbstractMatrix{T}
    sub_diagonal::Vector{SparseMatrixCSC{T, Int64}}
    main_diagonal::Vector{SparseMatrixCSC{T, Int64}}
    super_diagonal::Vector{SparseMatrixCSC{T, Int64}}
end

# Define an outer constructor which can infer the type
# required for creation of this matrix
function BlockTridiagonal(sub::Vector{SparseMatrixCSC{T, Int64}},
    main::Vector{SparseMatrixCSC{T, Int64}},
    super::Vector{SparseMatrixCSC{T, Int64}}) where T
return BlockTridiagonal{T}(sub, main, super)
end

# Define the size function for BlockTridiagonal to allow matrix operations
Base.size(B::BlockTridiagonal) = (length(B.main_diagonal) * size(B.main_diagonal[1], 1),
                                  length(B.main_diagonal) * size(B.main_diagonal[1], 2))


# Method to access the elements of the BlockTridiagonal matrix
# row, col = i,j
function Base.getindex(B::BlockTridiagonal{T}, i::Int, j::Int) where T
    block_size = size(B.main_diagonal[1], 1)
    block_row = (i + block_size - 1) ÷ block_size
    block_col = (j + block_size - 1) ÷ block_size

    local_i = i - (block_row - 1) * block_size
    local_j = j - (block_col - 1) * block_size

    if block_row == block_col
        return B.main_diagonal[block_row][local_i, local_j]
    elseif block_row == block_col - 1
        return B.super_diagonal[block_row][local_i, local_j]
    elseif block_row - 1 == block_col
        return B.sub_diagonal[block_col][local_i, local_j]
    else
        return zero(T)
    end
end                                  


# Overload the * operator for BlockTridiagonal matrices and vectors
# UNTESTED
function Base.:*(B::BlockTridiagonal, x::AbstractVector)
    # Check dimensions
    n_blocks = length(B.main_diagonal)
    block_size = size(B.main_diagonal[1], 1)
    if size(x, 1) != n_blocks * block_size
        throw(DimensionMismatch("Matrix and vector sizes do not match"))
    end

    # Perform the multiplication
    result = zeros(eltype(x), size(x, 1))
    for i in 1:n_blocks
        # Multiply main diagonal block
        result[(i-1)*block_size+1:i*block_size] += B.main_diagonal[i] * x[(i-1)*block_size+1:i*block_size]
        
        # Multiply super diagonal block if it exists
        if i < n_blocks
            result[(i-1)*block_size+1:i*block_size] += B.super_diagonal[i] * x[i*block_size+1:(i+1)*block_size]
        end
        
        # Multiply sub diagonal block if it exists
        if i > 1
            result[(i-1)*block_size+1:i*block_size] += B.sub_diagonal[i-1] * x[(i-2)*block_size+1:(i-1)*block_size]
        end 
    end
    
    return result
end

# Assuming `BlockTridiagonal` is already defined and `BT` is an instance of `BlockTridiagonal`
# UNTESTED
# unsure if needed
function block_tridiagonal_to_sparse(BT::BlockTridiagonal)
    # Combine the blocks from the BlockTridiagonal into a single SparseMatrixCSC
    # This is a simplified example and assumes that the blocks are correctly aligned and sized
    n_blocks = length(BT.main_diagonal)
    block_size = size(BT.main_diagonal[1], 1)
    full_size = n_blocks * block_size
    A = spzeros(full_size, full_size)

    for i in 1:n_blocks
        row_offset = (i-1) * block_size
        col_offset = (i-1) * block_size

        # Place main diagonal blocks
        A[row_offset+1:row_offset+block_size, col_offset+1:col_offset+block_size] = BT.main_diagonal[i]

        # Place super diagonal blocks
        if i < n_blocks
            A[row_offset+1:row_offset+block_size, col_offset+block_size+1:col_offset+2*block_size] = BT.super_diagonal[i]
        end

        # Place sub diagonal blocks
        if i > 1
            A[row_offset-block_size+1:row_offset, col_offset+1:col_offset+block_size] = BT.sub_diagonal[i-1]
        end
    end

    return A
end