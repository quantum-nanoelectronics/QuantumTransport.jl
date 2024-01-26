using LinearAlgebra
using SparseSuite

⊗(A,B) = kron(A,B)


σ₂ = [0 -im; im 0]

# here are some matrices you can use!
# makes a (n x m) x (n x m) matrix, where m is the block size.  
function effectiveMass(n::Int, m::Int=2)
        return Tridiagonal(ones(n-1),-(2+im*10^-5)*ones(n),ones(n-1))⊗I(m)
end

# block size 2
function spinOrbitHamiltonian(n::Int)
        return Tridiagonal(ones(n-1), 0*ones(n), ones(n-1))⊗σ₂ + -1*ones(n*2)
end

# this takes in your approximation Gʳ(Energy::Float64), the full LinAlg inv() Gʳ, and the size of the matrix. 
function verifyCorrectness(approximatedGʳ::Function, fullGʳ::Function, N::Int, Emin::Float64=-3.0, Emax::Float64=3.0, nE::Int=50, cutoff=0.05)
        sumdiff = 0
        Evals = LinRange(Emin,Emax, nE)
        for E ∈ Evals
               sumdiff += sum(abs.(abs.(approximatedGr(E)).*(approximatedGr(E) - fullGr(E))))/(n)^2 
        end
        sumdiff = sumdiff/nE
        if(sumdiff < cutoff)
                return true
        else
                return false
        end
end

function 
