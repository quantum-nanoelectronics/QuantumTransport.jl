using LinearAlgebra

"""
    generate_matrix()

Generate a matrix based on some criteria.

"""
function generate_matrix()
    A = rand(0:100, 1100, 1100)

    dimensions = size(A)
    if dimensions[1] != dimensions[2]
        error("Dimensions mismatch...")
    end

    if LinearAlgebra.det(A) == 0
        error("Det = 0, exiting...")
    end

    return A
end


"""
    block_inversion(matrix, threshold=8)

Compute the inverse of a block matrix using a recursive algorithm.

# Arguments
- `matrix`: The block matrix to be inverted.
- `threshold`: The maximum size of the blocks at which the algorithm switches to directly computing the inverse using `inv()`.

# Returns
The inverse of the input block matrix.

"""
function block_inversion(matrix, threshold=8)
    rows, cols = size(matrix)
    
    if max(rows, cols) <= threshold
        return inv(matrix)  # It will throw an error if the matrix is singular
    end
    
    half_rows, half_cols = div(rows, 2), div(cols, 2)
    A = matrix[1:half_rows, 1:half_cols]
    B = matrix[1:half_rows, half_cols+1:end]
    C = matrix[half_rows+1:end, 1:half_cols]
    D = matrix[half_rows+1:end, half_cols+1:end]
    
    A_inv = block_inversion(A, threshold)
    D_CA_inv_B = D - C * A_inv * B
    
    D_CA_inv_B_inv = block_inversion(D_CA_inv_B, threshold)
    
    upper_left = A_inv + A_inv * B * D_CA_inv_B_inv * C * A_inv
    upper_right = -A_inv * B * D_CA_inv_B_inv
    lower_left = -D_CA_inv_B_inv * C * A_inv
    lower_right = D_CA_inv_B_inv
    
    # Combining the results into a single matrix
    result = [upper_left upper_right; lower_left lower_right]
    return result
end


