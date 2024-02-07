#using Pkg; Pkg.add("Arpack");
#using Pkg; Pkg.add(Pkg.PackageSpec(;name="Arpack", version="0.5"));
using LinearAlgebra
using SparseArrays
#using Arpack: eigs
include("./Decomp.jl")

offset(A) = A .+ 1e-5im

function construct_laplace(n)
    #Construct a sparse diagonal matrix from Pairs of vectors and diagonals.
    #Each vector (second) will be placed on the (first) diagonal.
    return spdiagm(-1 => offset(ones(ComplexF64, n-1)),
                    0 => offset(-2*ones(ComplexF64,n)),
                    1 => offset(ones(ComplexF64,n-1)))
end

n = 4
H = construct_laplace(n)
test_struct = EigDecomp_mod.EigDecomp(Matrix(H))
