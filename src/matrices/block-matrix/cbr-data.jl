using LinearAlgebra
using SparseArrays
using Arpack

include("./block-tridiagonal.jl")

mutable struct Eigens{T}
    values::Vector{T}
    vectors::Matrix{T}
end


# Use ARPACK to compute the eigenvalues/vectors closest to the given energy level
function compute_closest_eigenvalues(BT::BlockTridiagonal, energy_level::Real, num_eigenvalues::Int)
    # Convert BlockTridiagonal to SparseMatrixCSC
    sparse_BT = block_tridiagonal_to_sparse(BT)
    # println(sparse_BT)

    # Use ARPACK's eigs function to find eigenvalues and eigenvectors
    # Set which=:LM because we are using shift-invert mode with a non-zero sigma
    return eigs(sparse_BT, nev=num_eigenvalues, which=:LM, sigma=energy_level)
end

# not needed
function print_blocks(label::String, blocks::Vector{<:SparseMatrixCSC})
    println("$label blocks:")
    for (i, block) in enumerate(blocks)
        println("Block $i:")
        display(block) # 'display' is used here for a better formatting of the sparse matrix
    end
end

# Create some sparse blocks here
block_size = 2
n_blocks = 3
A = [sprand(block_size, block_size, 0.5) for _ in 1:n_blocks]
B = [sprand(block_size, block_size, 0.5) for _ in 1:n_blocks-1]
C = [sprand(block_size, block_size, 0.5) for _ in 1:n_blocks-1]

# Example usage for A, B, and C:
print_blocks("A", A)
print_blocks("B", B)
print_blocks("C", C)
println(typeof(A)) # Vector{SparseMatrixCSC{Float64, Int64}}

# Create BlockTridiagonal matrices
BT = BlockTridiagonal(B, A, C)
BT1 = BlockTridiagonal(C, A, B)

# Testing multiplication by a vector
x = rand(n_blocks * block_size)
result = @which BT * x
println(result)

getindex(BT, 1, 7)
energy_level = 1e-5  # Interested energy level
num_eigenvalues = 1  # Number of eigenvalues/vectors

eigenvalues, eigenvectors = compute_closest_eigenvalues(BT, energy_level, num_eigenvalues)
e1 = Eigens{}(eigenvalues, eigenvectors)
println(e1)

