module LowRankMatrices
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod
    using Arpack

    export LRA
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
                Emin::𝐑, Emax::𝐑; 
                initialrank::Int = 64, rankstep::Int = 32,
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

            λs = undef; Rvecs = undef; Lvecs = undef
            curr_rank = initialrank
            success = false
            iters = 0
            centralenergy = (Emax+Emin)/2
            shiftedA = A - centralenergy*I(A.n)
            energywindow = Emax-Emin
            #= 
             Two stop conditions: number of iterations exceeds allowed maxiumum or
             eigen-energies are sufficiently close to defined window and the normalized
             error is small. In case maximum iterations achieved, will return original matrix.
            =#

            #= NOTE: second condition unnecessary, but useful for debugging small matrices =#
            while iters <= maxiters && curr_rank <= A.n
                # this needs to be changed, SR doesn't work -- needs SM, which is not implemented in arnoldi.jl
                # it is implemented in ARPACK.jl though
                #=schur, hist = partialschur(shiftedA, nev = curr_rank, tol = sqrt(eps()), which = ArnoldiMethod.SR())
                λs, Rvecs = partialeigen(schur) =#
                error = 0
                #error = norm_error(A, reconstruct(Rvecs, conj(Rvecs), λs)) #NOTE: expensive
                #println(error)
                λs, Rvecs = eigs(shiftedA, nev = curr_rank, tol = sqrt(eps()), which=:SM)
                _, Lvecs = eigs(shiftedA', nev = curr_rank, tol = sqrt(eps()), which=:SM)
                minλ = minimum(real.(λs))
                maxλ = maximum(real.(λs))
                println("Rank = $curr_rank, Min λ = $minλ, Max λ = $maxλ, Energy window = $energywindow")
                if minλ < - energywindow/2 && maxλ > energywindow/2
                    success = true
                    break
                end
                #=if abs(Emin - minλ) < ΔE &&
                   abs(Emax - maxλ) < ΔE &&
                   error < errortol
                    success = true
                    break
                end=#
                iters+=1; curr_rank+=rankstep
            end

            #if success return new(Rvecs, conj(Rvecs), λs, curr_rank, A.n) else return A end
            if success return new(Rvecs, Lvecs, λs.+centralenergy, curr_rank, A.n) else return nothing end

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
        print(io, "Left eigenvectors:\n")
        display(A′.Lvecs)
        #NOTE: can display if debugging
        print(io, "Reconstructed ∑λvvᴴ:\n")
        display(reconstruct(A′))
    end

    # world ending critical failure if reconstruct. reconstruct expects normal matrix
    function reconstruct(A′::LRA)
        A = zeros(ComplexF64, A′.n, A′.n)
        for k in eachindex(A′.λs)
            A .+= A′.λs[k] * A′.Rvecs[:,k] * transpose( A′.Lvecs[:,k] )
        end
        return A
    end

    function reconstruct(Rvecs::Matrix{ComplexF64}, Lvecs::Union{Matrix{ComplexF64}, Transpose{ComplexF64}},
                         λs::Array{ComplexF64})
        A = zeros(ComplexF64, size(Rvecs)[1], size(Rvecs)[1])
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
        Ω = zeros(ComplexF64, B′.n, B′.rank)   
        ω = Matrix{ComplexF64}(undef, B′.n, B′.rank)  
        λs = Matrix{ComplexF64}(undef, B′.rank, 1)   
        for j in eachindex(B′.λs)
            for i in eachindex(A′.λs)
                Ω[:,j] += (A′.λs[i] * B′.λs[j] * A′.Lvecs[:,i]' * B′.Rvecs[:,j]) * A′.Rvecs[:,i] 
            end
            #n = dot(conj(transpose(Ω[:,j])), Ω[:,j])
            n = Ω[:,j]' * Ω[:,j]
            λs[j] = n
            ω[:,j] = Ω[:,j] / n
        end
        C′ = LRA(ω, B′.Lvecs, λs, B′.rank, B′.n)
        return C′
    end


    function LRAbyMatrix(A′::LRA, A::Union{SparseMatrixCSC{ComplexF64, Int64}, Matrix{ComplexF64}})
        if A′.n != size(A)[2]
            throw(DomainError(A′, "BAD SHAPE IN LRA MATRIX MULTIPLICATION"))
        end
        #NOTE: Assumes A is linear
        C′ = LRA( A * A′.Rvecs, A′.Lvecs, A′.λs, A′.rank, A′.n)
        return C′
    end

    #TODO: maybe a nice operator representation of this?
    function mul_ind_LRA(A′::LRA, B′::LRA, i::Int, j::Int)
        return A′[i, :] ⋅ B′[:,j]
    end

    function add_scaled_identity(D::Diagonal, A′::LRA)
        D_1 = D[1]
        for i in eachindex(D)
            if D[i] != D_1
                throw(DomainError(A′, "Perturbative sum not yet implemented for Diagonal matrices"))
            end
        end
        A′.λs += D_1
    end

    function inv(A::LRA)
        invλs = (A.λs).^(-1)
        return LRA(A.Rvecs, A.Lvecs, invλs, A.rank, A.n)
    end

    function perturbative_sum(A::LRA, B::LRA)
        # figure out which matrix is bigger
        sumλs_A = sum(abs.(A.λs))
        sumλs_B = sum(abs.(B.λs))
        H = undef; δH = undef;
        if sumλs_A > sumλs_B
            H = A
            δH = B
        else
            δH = A
            H = B
        end
        # now do non-hermitian degenerate perturbation theory
        Heff = zeros(ComplexF64,H.rank,H.rank)
        Heff += Diagonal(H.λs)
        for row ∈ 1:H.rank
            for col ∈ 1:H.rank
                Heff[row, col] = Matrix(H.Lvecs[:,row]')*δH*H.Rvecs[:,col]
            end
        end
        # Rcoeffs and Lcoeffs are rank × n matrices
        λs, Rcoeffs = eigen(Heff)
        Lcoeffs = eigvecs(Heff')
        newRvecs = H.Rvecs*Rcoeffs
        newLvecs = H.Lvecs*Lcoeffs
        return LRA(newRvecs, newLvecs, λs, H.rank, H.n)
    end

    function number_multiplication(A′::LRA, C::Number)
        return LRA(A′.Rvecs, A′.Lvecs, C*A′.λs, A′.rank, A′.n)
    end

    function LRA_vector_multiplication(A::LRA, V::Vector)
        # assert size of vector same as A
        return A.Lvecs*(Diagonal(A.λs)*(A.Rvecs*V))
    end
 
    function LRA_vector_multiplication(V::Vector, A::LRA)
        # assert size of vector same as A
        return ((V*A.Lvecs)*Diagonal(A.λs))*A.Rvecs
    end
    
    Base.:size(A′::LRA) = A′.n, A′.n
    Base.:getindex(A′::LRA, i::Int, j::Int) = getindex_LRA(A′::LRA, i::Int, j::Int)
    
    #TODO: define rule for conversion?
    #Base.:promote_rule(...)
    #Base.:convert(...)

    # TODO: define shifting of eigenenergies when add identity matrix
    # TODO: define addition

    #Base.:+(A′::LRA, A::T) where T<:Matrix = +(promote(A′,A)...)
    #Base.:-(A′::LRA, A::T) where T<:Matrix = -(promote(A′,A)...)
    #Base.:+(A′::LRA, B′::LRA) = +(reconstruct(A′), reconstruct(B′))
    #Base.:-(A′::LRA, B′::LRA) = -(reconstruct(A′), reconstruct(B′))
    Base.:+(D::Diagonal, A′::LRA) = add_scaled_identity(D,A′)
    Base.:+(A′::LRA, D::Diagonal) = add_scaled_identity(D,A′)
    Base.:+(A::LRA,B::LRA) = perturbative_sum(A,B)
    Base.:inv(A::LRA) = inv(A)

    Base.:*(V::Vector,A′::LRA) = LRA_vector_multiplication(V, A′)  
    Base.:*(A′::LRA,V::Vector) = LRA_vector_multiplication(A′, V)
    Base.:*(A′::LRA,C::Number) = number_multiplication(A′, C)
    Base.:*(C::Number, A′::LRA) = number_multiplication(A′, C)
    Base.:*(A′::LRA, A::Union{SparseMatrixCSC{ComplexF64, Int64}, Matrix{ComplexF64}}) = LRAbyMatrix(A′, A)
    Base.:*(A::Union{SparseMatrixCSC{ComplexF64, Int64}, Matrix{ComplexF64}}, A′::LRA) = LRAbyMatrix(A′, A)
    Base.:*(A′::LRA, B′::LRA) = LRAbyLRA(A′, B′)
end
