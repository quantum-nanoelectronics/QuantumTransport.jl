using LinearAlgebra
using SparseArrays

⊗(A,B) = kron(A,B)
σ₂ = [0 -im; im 0]

# makes a (n x m) x (n x m) matrix, where m is the block size.
# this function was changed such that all elements are complex.
function effectiveMass(n::Int, m::Int=2)
        off_diag = ComplexF64.(ones(n-1))  # Converted off-diagonal elements to ComplexF64
        diag = -(2+im*10^-5)*ones(ComplexF64, n)  # Ensured diagonal elements are ComplexF64
        return Tridiagonal(off_diag, diag, off_diag) ⊗ I(m)
        # return Tridiagonal(ones(n-1),-(2+im*10^-5)*ones(n),ones(n-1))⊗I(m)
end

# Fixed this function by adding Tridiagonals of the same type
# Creates a tridiagonal with a block size of 2. 
function spinOrbitHamiltonian(n::Int)
        return Tridiagonal(ones(n-1), zeros(n), ones(n-1))⊗σ₂ + Tridiagonal(ComplexF64[(-1)*mod(i, 2) for i in 1:(2n-1)], ComplexF64[-1 for i in 1:(2n)], ComplexF64[(-1)*mod(i, 2) for i in 1:(2n-1)])
end


# do not use this for benchmarking
# creates a matrix with random hoppings
function randomHoppingHamiltonian(n::Int,m::Int=4)
        A = zeros(ComplexF64,n*m,n*m)
        A += 1E-4*im*I(m*n)
        for i = 1:(n-1)
                colvec = zeros(n)
                rowvec = zeros(n)'
                rowvec[i] = 1;
                colvec[i+1] = 1;
                t = rand(ComplexF64,m,m)
                # do the top 
                offdiagBlock = colvec*rowvec⊗t
                A += offdiagBlock + offdiagBlock'
        end
        return A
end



function verifyCorrectness(approximatedGʳ::Function, fullGʳ::Function, N::Int, Emin::Float64=-3.0, Emax::Float64=3.0, nE::Int=50, cutoff=0.05)
        sumdiff = 0
        Evals = LinRange(Emin,Emax, nE)
        for E ∈ Evals
               sumdiff += sum(abs.(abs.(approximatedGʳ(E)).*(approximatedGʳ(E) - fullGʳ(E))))/(N)^2 
        end
        sumdiff = sumdiff/nE
        println("Error = ", sumdiff)
        if(sumdiff < cutoff)
                return true
        else
                return false
        end
end