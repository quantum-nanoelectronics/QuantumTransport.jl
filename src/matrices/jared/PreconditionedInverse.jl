module PreconditionedInverse
    using LinearAlgebra
    using SparseArrays
    using ArnoldiMethod
    using Arpack
    include("LowRankMatrices.jl")
    using .LowRankMatrices
    
    export eigen

    function eigen(matrixsum::Tuple{Vararg}, 
        guess_eigvecs::Union{Matrix{ComplexF64},Matrix{Float64}}, 
        guess_eigvals::Union{Vector{Float64},Vector{ComplexF64}}, tolerance::Float64=1E-5)
        
        rank = size(guess_eigvals)[1]; n = size(guess_eigvecs)[1]
        final_Reigvecs = zeros(ComplexF64,n,rank)
        final_Leigvecs = zeros(ComplexF64,rank,n)
        final_Reigvals = zeros(ComplexF64,rank)
        final_Leigvals = zeros(ComplexF64,rank)
        for i = 1:rank
            guess_Reigvec = guess_eigvecs[:,i]
            guess_Leigvec = guess_eigvecs[:,i]'
            iter_Reigvec = zeros(ComplexF64,n,1)
            iter_Leigvec = zeros(ComplexF64,n,1)
            error = 1
            while error > tolerance
                for M ∈ matrixsum
                    iter_Reigvec += M*guess_Reigvec
                    iter_Leigvec += guess_Leigvec*M
                end
                guess_Reigvec = iter_Reigvec/norm(iter_Reigvec)
                guess_Leigvec = iter_Leigvec/norm(iter_Leigvec)
                error = (guess_Reigvec - iter_Reigvec)'*(guess_Reigvec -  iter_Reigvec) + 
                (guess_Leigvec - iter_Leigvec)*(guess_Leigvec -  iter_Leigvec)'
                println("error = $error")
            end
            λR = norm(iter_Reigvec)
            λL = norm(iter_Leigvec)
            final_Reigvec = iter_Reigvec / norm(iter_Reigvec)
            final_Leigvec = iter_Leigvec / norm(iter_Leigvec)
            final_Reigvec = final_Reigvec/(final_Leigvec*final_Reigvec)
            final_Reigvecs[:,i] = final_Reigvec
            final_Leigvecs[i,:] = final_Leigvec
            final_Reigval[i] = λR
            final_Leigval[i] = λL
        end
        return final_Leigvecs, final_Leigvals, final_Reigvals, final_Reigvecs
    end
    #Base.:eigen(matrixsum::Tuple{AbstractMatrix}, guess_eigvecs::Union{Matrix{ComplexF64},Matrix{Float64}}, guess_eigvals::Union{Vector{Float64},Vector{ComplexF64}}, tolerance::Float64=1E-5) = eigen(matrixsum, guess_eigvecs, guess_eigvals, tolerance)
end