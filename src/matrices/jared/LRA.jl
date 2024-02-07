#using Pkg; Pkg.add("Arpack");
using LinearAlgebra
using SparseArrays
#using Arpack: eigs
using ArnoldiMethod


offset(A) = A .+ 1e-5im

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

# takes nxn size
# returns ComplexF64 laplace operator
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

module eigen_decompose
    using LinearAlgebra
    using SparseArrays
    using Arpack: eigs
    #using ArnoldiMethod

    export EigenDecomp

    struct EigenDecomp #<: AbstractArray{XXX, 2}
        #F :: Eigen 
        e_vals :: ComplexF64
        e_vecs :: Array
        function EigenDecomp(A)#(M:Array{XXX, 2})
            #decomp, history = partialschur(A, nev=10, tol=1e-6, which=SR());
            #e_vals, e_vecs = partialeigen(decomp)
            e_vals, e_vecs = eigs(A, nev=10)
            #e_vecs = eigvecs(A)
        end
    end

    function Base.show(io::IO, M::EigenDecomp)
        # type inference cannot occur here
        print(io,M.e_vals,"\n")
        print(io,M.e_vecs,"\n")
    end
end

n = 14
H = construct_laplace(n)

test_struct = eigen_decompose.EigenDecomp(Matrix(H))

"
#display(H)
# currently converting from sparse 
E = eigvals(Tridiagonal(H))
ϕ = eigvecs(Tridiagonal(H))
#eigs(H, nev = 2, which=:SM)
#display(ϕ)
E, ϕ = truncate_eigens(E, ϕ, 3)
S = deepcopy(ϕ)
X = construct_pos(n)
Sꜛ = conj(transpose(S))
Xₛ = Sꜛ * X * S
xₛ = eigvals(Tridiagonal(Xₛ))
ψ  = eigvecs(Tridiagonal(Xₛ))
#Iₙ = I(n)
E_iη = offset(E)
Gᵣ = inv( Matrix(Diagonal(E_iη) - H) )
"
