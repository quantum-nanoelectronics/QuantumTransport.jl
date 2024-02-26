using LinearAlgebra
using SparseArrays
using ArnoldiMethod
include("./LRA_type.jl")

offset(A) = A .+ 1e-5im

function construct_laplace(n)
    #Construct a sparse diagonal matrix from Pairs of vectors and diagonals.
    #Each vector (second) will be placed on the (first) diagonal.
    return spdiagm(-1 => offset(ones(ComplexF64, n-1)),
                    0 => offset(-2*ones(ComplexF64,n)),
                    1 => offset(ones(ComplexF64,n-1)))
end

# takes nxn size
# returns ComplexF64 position operator
function construct_pos(n)
    return spdiagm(-1 => (sqrt.(1:n-1)),
                    1 => (sqrt.(1:n-1)))
end

function truncate_eigens(eigen_vals, eigen_vecs, cutoff)
    for i in eachindex(eigen_vals)
        # complex magnitude abs()
        if abs(eigen_vals[i]) > cutoff
            eigen_vecs[:, i] .= 0
            eigen_vals[i] = 0
        end
    end
    return eigen_vals, eigen_vecs
end


X = SparseMatrixCSC(  [ 1.0 1.0im; -1.0im 1.0 ] )
display(X)
L = LRA_mod.LRA(X, 0.3, 3.0, 0.25)
print(L)

