module LRA_mod
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod

    # can likely pull these functions in from elsewhere
    norm_error(A,B) = norm( A - B ) / norm(A)
    
    #= in -> sparse, square matrix of value type ComplexF64 =#
    struct LRA #<: AbstractMatrix{ComplexF64}
        Rvecs :: Matrix{ComplexF64}
        Lvecs :: Matrix{ComplexF64}
        λs    :: Array{ComplexF64}
        rank  :: Int64
        n     :: Int64

        function LRA(A::SparseMatrixCSC{ComplexF64, Int64},
                Emin::𝐑, Emax::𝐑, ΔE::𝐑,
                initialrank::Int = 10, rankstep::Int = 1,
                maxiters::Int = 100, errortol::𝐑 = 1e-2) where 𝐑 <: Real
            if A.n != A.m
                throw(DomainError(A, "LRA DEFINED FOR SQUARE MATRICES ONLY."))
            end

            #= Only normal matrices are supported s.t. vₗ = vᵣ* =#
            #= TODO: For hermitian matrices, implement Lanczos method =#
             if norm_error(A * conj(A) , conj(A) * A) > errortol
                throw(DomainError(A, "LRA NOT DEFINED FOR NON NORMAL MATRIX YET"))
            end

            #= Ignore small matrices. Decrease INITIAL_RANK for testing. =#
            if A.n < initialrank
                return A
            end

            λs = undef; Rvecs = undef
            curr_rank = initialrank
            success = false
            iters = 0

            #= 
             Two stop conditions: number of iterations exceeds allowed maxiumum or
             eigen-energies are sufficiently close to defined window and the normalized
             error is small. In case maximum iterations achieved, will return original matrix.
            =#

            #= NOTE: second condition unnecessary, but useful for debugging small matrices =#
            while iters <= maxiters && curr_rank <= A.n
                schur, hist = partialschur(A, nev = curr_rank, tol = sqrt(eps()), which = ArnoldiMethod.SR())
                λs, Rvecs = partialeigen(schur)
                error = norm_error(A, reconstruct(Rvecs, conj(Rvecs), λs)) #NOTE: expensive
                if abs(Emin - minimum(real.(λs))) < ΔE &&
                   abs(Emax - maximum(real.(λs))) < ΔE &&
                   error < errortol
                    success = true
                    break
                end
                iters+=1; curr_rank+=rankstep
            end

            if success return new(Rvecs, conj(Rvecs), λs, curr_rank, A.n) else return A end

        end

        function LRA(Rvecs::Matrix{ComplexF64}, Lvecs::Matrix{ComplexF64},
                     λs::Array{ComplexF64}, rank::Int64, n::Int64)
            return new(Rvecs, Lvecs, λs, rank, n)
        end

    end

    function Base.show(io::IO, A′::LRA)
        print(io, "Decomposed matrix of rank ", A′.rank, "\n")
        print(io, "Eigenvalues:\n")
        display(A′.λs)
        print(io, "Right eigenvectors:\n")
        display(A′.Rvecs)
        #NOTE: can display if debugging
        print(io, "Reconstructed ∑λvvᴴ:\n")
        display(reconstruct(A′))
    end

    # world ending critical failure if reconstruct. reconstruct expects normal matrix
    function reconstruct(A′::LRA)
        A = Matrix{ComplexF64}(undef, A′.n, A′.n)
        for k in eachindex(A′.λs)
            A .+= A′.λs[k] * A′.Rvecs[:,k] * transpose( A′.Lvecs[:,k] )
        end
        return A
    end

    function reconstruct(Rvecs::Matrix{ComplexF64}, Lvecs::Union{Matrix{ComplexF64}, Transpose{ComplexF64}},
                         λs::Array{ComplexF64})
        A = Matrix{ComplexF64}(undef, size(Rvecs)[1], size(Rvecs)[1])
        for k in eachindex(λs)
            A .+= λs[k] * Rvecs[:,k] * transpose( Lvecs[:,k] )
        end
        return A
    end

    function getindex_LRA(A′::LRA, i::Int, j::Int)
        Aᵢⱼ = 0
        for k in eachindex(A′.λs)
            Aᵢⱼ += A′.λs[k] * A′.Rvecs[i,k] * A′.Lvecs[j,k]
        end
        return Aᵢⱼ
    end

    function trace(A′::LRA)
        return sum(A′.λs)
    end

    function LRAbyLRA(A′::LRA, B′::LRA)
        if A′.n != B′.n
            throw(DomainError(A′, "BAD SHAPE IN LRA MATRIX MULTIPLICATION"))
        end
        #FIXME: still a matrix product
        #NOTE: Assumes A is linear
        λs = ( A′.Rvecs ⋅ B′.Rvecs ) * ( A′.λs ⋅ B′.λs )
        C′ = LRA( reconstruct(A′) * B′.Rvecs, B′.Lvecs, B′.λs, B′.rank, B′.n)
        return C′
    end

    #=
    function LRAbyLRA(A′::LRA, B′::LRA)
        if A′.n != B′.n
            throw(DomainError(A′, "BAD SHAPE IN LRA MATRIX MULTIPLICATION"))
        end
        #FIXME: still a matrix product
        #NOTE: Assumes A is linear
        C′ = LRA( reconstruct(A′) * B′.Rvecs, B′.Lvecs, B′.λs, B′.rank, B′.n)
        return C′
    end
    =#

    function LRAbyMatrix(A′::LRA, A::Union{SparseMatrixCSC{ComplexF64, Int64}, Matrix{ComplexF64}})
        if A′.n != A.n
            throw(DomainError(A′, "BAD SHAPE IN LRA MATRIX MULTIPLICATION"))
        end
        C′ = LRA( A * A′.Rvecs, A′.Lvecs, A′.λs, A′.rank, A′.n)
        return C′
    end

    #TODO: maybe a nice operator representation of this?
    function mul_ind_LRA(A′::LRA, B′::LRA, i::Int, j::Int)
        return A′[i, :] ⋅ B′[:,j]
    end

    Base.:size(A′::LRA) = A′.n, A′.n
    Base.:getindex(A′::LRA, i::Int, j::Int) = getindex_LRA(A′::LRA, i::Int, j::Int)
    
    #TODO: convert down to LRA?
    #Base.:promote_rule(...)
    #Base.:convert(...)

    # TODO: define shifting of eigenenergies when add identity matrix
    #Base.:+(A′::LRA, A::T) where T<:Matrix = +(promote(A′,A)...)
    #Base.:-(A′::LRA, A::T) where T<:Matrix = -(promote(A′,A)...)
    #Base.:+(A′::LRA, B′::LRA) = +(reconstruct(A′), reconstruct(B′))
    #Base.:-(A′::LRA, B′::LRA) = -(reconstruct(A′), reconstruct(B′))

    Base.:*(A′::LRA, A::Union{SparseMatrixCSC{ComplexF64, Int64}, Matrix{ComplexF64}}) = LRAbyMatrix(A′, A)
    Base.:*(A::Union{SparseMatrixCSC{ComplexF64, Int64}, Matrix{ComplexF64}}, A′::LRA) = LRAbyMatrix(A′, A)
    Base.:*(A′::LRA, B′::LRA) = LRAbyLRA(A′, B′)
end
