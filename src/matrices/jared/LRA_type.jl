module LRA_mod
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod
    # ArnoldiMethod implements eigs from ARPACK... ARPACK machine broke
    # NOTE: for hermitian matrices, the Arnoldi iteration reduces to the
    # Lanczos algorithm — potential speed improvement

    # in  -> sparse matrix
    # holds the deconstructed, low rank matrix
    struct LRA <: AbstractMatrix{ComplexF64}
        R :: Matrix{ComplexF64}
        #L :: Matrix{ComplexF64}
        λ :: Array{ComplexF64}
        Λ :: Matrix{ComplexF64}
        rank :: Int64
        shape :: Tuple{Int64, Int64}

        function LRA(A::SparseMatrixCSC{ComplexF64, Int64},
                           Emin::Rt, Emax::Rt, delta::Rt) where Rt <: Real
            if A.n != A.m
                throw(DomainError(A, "LRA DEFINED FOR SQUARE MATRICES ONLY."))
            end

            if A * conj(A) == conj(A) * A
                #sub hermition condiion for able to get L = R*
                last_rank = 0 
                current_rank = 0
                λ = undef
                R = undef
                while current_rank <= A.n
                    last_rank = Int(ceil(sqrt(A.n)))
                    current_rank = last_rank
                    R_decomp, R_hist = partialschur(A, nev = current_rank, tol = sqrt(eps()), which = ArnoldiMethod.SR())
                    λ, R = partialeigen(R_decomp)
                    if ((abs(Emin - minimum(real.(λ))) < delta && abs(Emax - maximum(real.(λ))) < delta)
                        || (minimum(real.(λ)) < Emin && maximum(real.(λ)) > Emax))
                       break
                    end
                    current_rank = Int(ceil(sqrt(A.n)))
                    current_rank = (current_rank <= A.n) ? current_rank*2 : A.n
                end
            else
                throw(DomainError(A, "LRA NOT DEFINED FOR NON NORMAL MATRIX YET"))
                #R_decomp, R_hist = partialschur(A, nev = A.n, tol = sqrt(eps()), which = SR())
                #L_decomp, L_hist = partialschur(transpose(A), nev = A.n, tol = sqrt(eps()), which = SR())
            end
            new(R, λ, diagm(λ), last_rank, size(A))
        end
    end

    function Base.show(io::IO, A::LRA)
        print(io, "Decomposed matrix of rank ", A.rank, "\n")
        print(io, "Eigenvalues:\n")
        display(A.λ)
        print(io, "Right eigenvectors:\n")
        display(A.R)
        print(io, "Reconstructed ∑λvvᴴ:\n")
        display(reconstruct(A))
    end

    # world ending critical failure if reconstruct
    function reconstruct(lra::LRA)
        A = Matrix{ComplexF64}(undef, lra.shape...)
        if A * conj(A) == conj(A) * A
            for j in eachindex(lra.λ)
                A .+= lra.λ[j] * lra.R[j,:] * transpose(conj(lra.R)[j,:])
            end
        else
            throw(DomainError(A, "LRA NOT DEFINED FOR NON NORMAL MATRIX YET"))
        end
        return A
    end

    #function reconstruct(lra::LRA)
    #    A = Matrix{ComplexF64}(undef, lra.shape...)
    #    if A * conj(A) == conj(A) * A
    #        return lra.R * lra.Λ * inv(lra.R)
    #    else
    #        throw(DomainError(A, "LRA NOT DEFINED FOR NON NORMAL MATRIX YET"))
    #    end
    #end

    #requires a better way to get index
    #function prod_element(A::T, B::T, i::I, j::I) where T<:LRA where I <: integer
    #    return A[i, :] ⋅ B[:,j]
    #end


    #Base.:getindex(A::T, i::I, j::I) where T<:LRA where I <: Integer = getindex_LRA(A, i, j)
    Base.:getindex(A::T) where T<:LRA = getindex(reconstruct(A))
    Base.:size(A::LRA) = A.shape
    Base.:promote_rule(::Type{Matrix{ComplexF64}}, ::Type{<:LRA}) = Matrix{ComplexF64}
    # convert up instead of down
    Base.:convert(::Type{<:Matrix{ComplexF64}}, A::LRA) = reconstruct(A)
    # define shifting of eigenenergies when add identity matrix

    Base.:+(x::LRA, y::T) where T<:Matrix = +(promote(x,y)...)
    Base.:*(x::LRA, y::T) where T<:Matrix = *(promote(x,y)...)
    Base.:-(x::LRA, y::T) where T<:Matrix = -(promote(x,y)...)
    Base.:/(x::LRA, y::T) where T<:Matrix = /(promote(x,y)...)
    Base.:+(x::LRA, y::LRA) = +(reconstruct(x), reconstruct(y))
    Base.:*(x::LRA, y::LRA) = *(reconstruct(x), reconstruct(y))
    Base.:-(x::LRA, y::LRA) = -(reconstruct(x), reconstruct(y))
    Base.:/(x::LRA, y::LRA) = /(reconstruct(x), reconstruct(y))
end
