using LinearAlgebra
using SparseArrays
using ArnoldiMethod
include("./LRA_type.jl")

#offset(A) = A .+ 1e-5im
offset(A) = A

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


#X = SparseMatrixCSC(  [ 1.0 1.0im; -1.0im 1.0 ] )
X = [ -1.0 0 0 -2im; 0 -1 2im 0; 0 -2im -1 0; 2im 0 0 -1 ]
Y = construct_laplace(4)
X = deepcopy(Y)

display(X)
display(Y)
#L = LRA_mod.LRA(X, -3.8, -0.38, 0.2, 2,1)
#L2 = deepcopy(L)
#X = LRA_mod.LRA(sparse(X), -2.9, 0.9, 1.0, 4, 1, 100, 1.0)
X = LRA_mod.LRA(Y, -3.8, -0.38, 0.25, 4, 1, 100, 1e-2)
Y = LRA_mod.LRA(Y, -3.8, -0.38, 0.25, 4, 1, 100, 1e-2)

Res = X * Y
println("calculated")
println(Res)

X_2 = construct_laplace(4) * construct_laplace(4)
X_2= LRA_mod.LRA(X_2, -0.77, 12.2, 1.0, 4, 1, 100, 1e-2)
println("correct===")
println(X_2)
