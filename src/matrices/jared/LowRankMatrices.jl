module LowRankMatrices
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod
    using Arpack

    export LRA
    # can likely pull these functions in from elsewhere
    norm_error(A,B) = norm( A - B ) / norm(A)
    
    #= in -> sparse, square matrix of value type ComplexF64 =#
    struct LRA <: AbstractMatrix{Complex}
        Rvecs :: Matrix{ComplexF64}
        Lvecs :: Matrix{ComplexF64}
        λs    :: Array{ComplexF64}
        rank  :: Int64
        n     :: Int64
        eigendecomposed :: Bool
        # only use this initializer with Hermitian sparse matrices
        function LRA(A::Union{SparseMatrixCSC{ComplexF64,Int64},SparseMatrixCSC{Float64,Int64}},
                Emin::𝐑, Emax::𝐑; 
                initialrank::Int = 64, rankstep::Int = 32,
                maxiters::Int = 20, errortol::𝐑 = 1e-2) where 𝐑 <: Real
            if A.n != A.m
                throw(DomainError(A, "LRA DEFINED FOR SQUARE MATRICES ONLY."))
            end
            #= Only normal matrices are supported s.t. vₗ = vᵣ* =#
            # TODO: add in a check for hermitian matrices specifically
            #= TODO: For hermitian matrices, implement Lanczos method =#
            #=
             if norm_error(A * conj(A) , conj(A) * A) > errortol
                throw(DomainError(A, "LRA NOT DEFINED FOR NON NORMAL MATRIX YET"))
            end
            =# 
            #= Ignore small matrices. Decrease INITIAL_RANK for testing. =#
            #=
            if A.n < initialrank
                return A
            end
            =#
            λs = undef; Rvecs = undef;
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
                λs, Rvecs = eigs(shiftedA, maxiter=1000, nev = curr_rank, which=:SM)
                #λs, Rvecs = eigs(shiftedA, maxiter=1000, nev = curr_rank, tol = sqrt(eps()), which=:SM)
                minλ = minimum(real.(λs))
                maxλ = maximum(real.(λs))
                # increase ARPACK decreases # eigvals
                curr_rank = length(λs)
                println("Rank = $curr_rank, Min λ = $minλ, Max λ = $maxλ, Energy window = $energywindow")
                if minλ < - energywindow/2 || maxλ > energywindow/2 
                    success = true
                    break
                end
                # normalize the phase on the first index of each eigvec
                #=if abs(Emin - minλ) < ΔE &&
                   abs(Emax - maxλ) < ΔE &&
                   error < errortol
                    success = true
                    break
                end=#
                iters+=1; curr_rank+=rankstep
            end
            #if success return new(Rvecs, conj(Rvecs), λs, curr_rank, A.n) else return A end
            return new(Rvecs, Array(Rvecs'), λs.+centralenergy, curr_rank, A.n, true)
        end

        function LRA(Rvecs::Matrix{ComplexF64}, Lvecs::Matrix{ComplexF64},
                     λs::Array{ComplexF64}, rank::Int64, n::Int64, eigendecomp::Bool=false)
            return new(Rvecs, Lvecs, λs, rank, n, eigendecomp)
        end

    end

    
    function Base.show(A::LRA)
        n = A.n
        println("$n×$n Low Rank Matrix of rank $(A.rank)")
        println("Eigenvalues: $(A.λs)")
        println("Size of left eigenvector block: $(size(A.Lvecs))")
        println("Size of right eigenvector block: $(size(A.Rvecs))")
    end

    function Base.show(io::IO, A::LRA)
        n = A.n
        println("$n×$n Low Rank Matrix of rank $(A.rank)")
        println("Eigenvalues: $(A.λs)")
        println("Size of left eigenvector block: $(size(A.Lvecs))")
        println("Size of right eigenvector block: $(size(A.Rvecs))")
    end

    # world ending critical failure if reconstruct. reconstruct expects normal matrix
    function reconstruct(A::LRA)
        @warn "You probably don't want to reconstruct a dense matrix. Reconsider, and do not use with large matrices."
        return A.Rvecs*Diagonal(A.λs)*A.Lvecs
    end

    function reconstruct(Rvecs::Matrix{ComplexF64}, Lvecs::Union{Matrix{ComplexF64}, Transpose{ComplexF64}},
                         λs::Array{ComplexF64})
        @warn "You probably don't want to do this. Reconsider, and do not use with large matrices."
        return Rvecs*Diagonal(A.λs)*Lvecs
    end

    function getindex_LRA(A::LRA, i::Int, j::Int)
        return A.Rvecs[i,:]*Diagonal(A.λs)*A.Lvecs[:,j]
    end

    function trace(A::LRA)
        sum = undef
        if A.eigendecomposed==true
            return sum(A.λs)
        else
            sum = 0 
            for i = 1:A.n
                sum += getindex_LRA(A,i,i)
            end
        end
        return sum
    end
    # TODO: fix, not a high priority though. 
    function LRAbyLRA(A::LRA, B::LRA)
        if A.n != B.n
            throw(DomainError(A, "BAD SHAPE IN LRA MATRIX MULTIPLICATION"))
        end
        # we wish to get RᴬΛᴬLᴬ*RᴮΛᴮLᴮ into the form RᶜΛᶜLᶜ = Rᴬ*RPᶜ*Λᶜ*LPᶜ*Lᴮ
        ΛᴬLᴬRᴮΛᴮ = Diagonal(A.λs)*A.Lvecs*B.Rvecs*Diagonal(B.λs)
        # this is done with inefficient, ill-conditioned jank because julia does not implement a left eigvec calculation. 
        F = eigen(ΛᴬLᴬRᴮΛᴮ); λsᶜ = F.values; RPᶜ = F.vectors
        LPᶜ = Base.inv(RPᶜ)   
        Rᶜ = A.Rvecs*RPᶜ; Lᶜ = LPᶜ*B.Lvecs
        return LRA(Rᶜ,Lᶜ,λsᶜ,B.rank, B.n, true)
    end
#=    function LRAbyLRA(A′::LRA, B′::LRA)
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
    end=#


    function LRAbyMatrix(A::LRA, B::M) where M <: AbstractMatrix
        if A.n != size(A)[2]
            throw(DomainError(A, "BAD SHAPE IN LRA MATRIX MULTIPLICATION"))
        end
        #NOTE: Assumes A is linear
        C = LRA(B*A.Rvecs, A.Lvecs, A.λs, A.rank, A.n, false)
        return C
    end

    #TODO: maybe a nice operator representation of this?
    function mul_ind_LRA(A′::LRA, B′::LRA, i::Int, j::Int)
        return A′[i, :] ⋅ B′[:,j]
    end

    function add_scaled_identity(D::Diagonal, A::LRA)
        D_1 = D[1,1]
        for i in 1:size(D)[1]
            if D[i,i] != D_1
                return perturbative_sum(A::LRA, D::Diagonal)
            end
        end
        return LRA(deepcopy(A.Rvecs), deepcopy(A.Lvecs), deepcopy(A.λs).+D_1, deepcopy(A.rank), deepcopy(A.n), deepcopy(A.eigendecomposed))
    end

    function inv(A::LRA)
        invλs = (A.λs).^(-1)
        return LRA(A.Rvecs, A.Lvecs, invλs, A.rank, A.n)
    end

    function perturbative_sum(A::LRA, B::M) where M <: AbstractMatrix 
        # now do non-hermitian degenerate perturbation theory
        return perturbative_sum(A,(B))
        return LRA(newRvecs, newLvecs, λs, A.rank, A.n, true)
    end


    function perturbative_sum(A::LRA, B::M...) where M <: AbstractMatrix 
        # now do non-hermitian degenerate perturbation theory
        Heff = zeros(ComplexF64,A.rank,A.rank)
        Heff += Diagonal(A.λs)
        for δH ∈ B 
            if typeof(δH) == LRA
                Test = A.Lvecs*δH*A.Rvecs
                Heff += A.Lvecs*δH.Rvecs*Diagonal(δH.λs)*δH.Lvecs*A.Rvecs
            else
                Test = A.Lvecs*δH*A.Rvecs
                Heff += Test
            end
        end
        # Rcoeffs and Lcoeffs are rank × n matrices

        F = eigen(Heff)
        λs = F.values; Rcoeffs = F.vectors
        Lcoeffs = Base.inv(Rcoeffs)
        newRvecs = A.Rvecs*Rcoeffs
        newLvecs = Lcoeffs*A.Lvecs
        return LRA(newRvecs, newLvecs, λs, A.rank, A.n, true)
    end

    function number_multiplication(A::LRA, C::Number)
        return LRA(A.Rvecs, A.Lvecs, C*A.λs, A.rank, A.n)
    end

    function LRA_vector_multiplication(A::LRA, V::Vector)
        # assert size of vector same as A
        return A.Rvecs*(Diagonal(A.λs)*(A.Lvecs*V))
    end
 
    function LRA_vector_multiplication(V::Vector, A::LRA)
        # assert size of vector same as A
        return ((V*A.Rvecs)*Diagonal(A.λs))*A.Lvecs
    end
    
    Base.:size(A′::LRA) = (A′.n, A′.n)
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
    Array(A::LRA) = reconstruct(A)
    Base.:+(D::Diagonal, A::LRA) = add_scaled_identity(D,A)
    Base.:+(A::LRA, B::AbstractMatrix...) = perturbative_sum(A, B)
    Base.:+(A::LRA, D::Diagonal) = add_scaled_identity(D,A)
    Base.:+(A::LRA,B::LRA) = perturbative_sum(A,B)
    Base.:inv(A::LRA) = inv(A)
    Base.display(A::LRA) = show(A)

    Base.:*(V::Vector,A′::LRA) = LRA_vector_multiplication(V, A′)  
    Base.:*(A′::LRA,V::Vector) = LRA_vector_multiplication(A′, V)
    Base.:*(A′::LRA,C::Number) = number_multiplication(A′, C)
    Base.:*(C::Number, A′::LRA) = number_multiplication(A′, C)
    Base.:*(A′::LRA, A::Union{SparseMatrixCSC{ComplexF64, Int64}, Matrix{ComplexF64}}) = LRAbyMatrix(A′, A)
    Base.:*(A::Union{SparseMatrixCSC{ComplexF64, Int64}, Matrix{ComplexF64}}, A′::LRA) = LRAbyMatrix(A′, A)
    Base.:*(A′::LRA, B′::LRA) = LRAbyLRA(A′, B′)
end