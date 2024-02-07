#+(x::M , y::ξ) = MyNumber(x.value + y.value)
# define all matrix operators on structure:
# +
# -
# kronecker
# cross
# dot
# ...

module EigDecomp_mod
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod
    #implements eigs from ARPACK... ARPACK machine broke
    #NOTE: for hermitian matrices, the Arnoldi iteration reduces to the
    # Lanczos algorithm — potential speed improvement

    #in-> sparse matrix
    #out'get' -> truncated sparse matrix from eigen decomp
    mutable struct EigDecomp{Tv<:Array{Complex}, Ti<:Complex} <: SparseArrays.AbstractSparseMatrix{Tv, Ti}
        v :: AbstractMatrix{Tv}
        λ :: AbstractArray{Tv,1}
        Q :: AbstractMatrix{Tv}
        Λ :: AbstractMatrix{Tv}
        history :: AbstractString

        function EigDecomp{Tv}(A::Tv) where Tv <: Array{Complex}
            decomp, history = partialschur(A, nev=10, tol=1e-6, which=SR());
            λ, v = partialeigen(decomp)
            Q = vcat[column(x) for x in v]
            Λ = Diagonal(λ)
        end
        #define constructor for non sparse matrix
        get_schur(A::EigDecomp) = return( A.Q*A.Λ*inv(A.Q) )
        get_herm(A::EigDecomp) = return( hcat( Λ[j] * v[j] * conj(transpose(v[j]))) )

        function Base.show(io::IO, M::EigDecomp)
            # type inference cannot occur here
            print(io, "Decomposed Matrix:\n")
            print(io, history, "\n")
            #print(io,get(M),"\n")
        end
    end
end
