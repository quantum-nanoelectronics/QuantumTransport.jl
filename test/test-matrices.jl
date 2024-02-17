# This file is included in matrix.jl

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

# Create/set test matrix
function setTestMatrix()
    # This is the initially created test matrix
    # return CreateBlockMatrix(args[1], args[2], args[3], args[5])  
    # effective mass test matrix.
    # return CreateBlockMatrix(args[1], args[2], args[3], args[5], sparse(effectiveMass(args[1] ÷ args[2], args[2])))
    # spin orbit hamiltonian test matrix - Block sizes of 2s only
    return CreateBlockMatrix(args[1], args[2], args[3], args[5], sparse(spinOrbitHamiltonian(args[1] ÷ args[2])))
end


function diagApproximatedGʳ(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + args[4]) * I - testMatrix.matrix
    return getInvRGFDiagonal(matrixCopy).matrix
end

function approximatedGʳ(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + args[4]) * I - testMatrix.matrix
    return getInvRGF(matrixCopy).matrix
end

function fullGʳ(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + args[4]) * I - testMatrix.matrix
    return getInvJulia(matrixCopy).matrix
end


# Function to time the two methods of computing the inverse of a block matrix.
function timeInv(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + args[4]) * I - testMatrix.matrix
    # Time the computation of the inverse of the dense matrix using built-in inversion
    print("Julia Inverse Time: ")
    @time juliaInv = getInvJulia(matrixCopy)
    # Time the computation of the inverse of the block matrix using RGF method
    print("RGF Inverse Time: ")
    @time rgfInv = getInvRGF(matrixCopy)
end






# other functions to test the correctness

function debugAllValues()
    a = fullGʳ(0.0)[:]
    b = approximatedGʳ(0.0)[:]
    for i in 1:size(a)[1]
            println(a[i], " ", b[i])
    end
end

function debugValues()
    # Print out comparisons between the computed blocks and the corresponding blocks in the dense inverse matrix
    # Comparing Diagonals
    println("Comparing Diagonals")
    # println(length(matrix[i:(i+numBlocks), i:(i+numBlocks)][:]))
    # println(length(diag[:]))
    for i in 1:matrixObject.numBlocks
            start_index = (i - 1) * matrixObject.blockSize + 1
            end_index = i * matrixObject.blockSize
            println(matrix[start_index:end_index, start_index:end_index])
            println(diag[i][:])
            println()
    end
    println()
    println()
    println()

    # Comparing Top Row
    println("Comparing Top")
    for i in 1:matrixObject.numBlocks
            start_index = (i - 1) * matrixObject.blockSize + 1
            end_index = i * matrixObject.blockSize
            println(matrix[1:matrixObject.blockSize, start_index:end_index])
            println(top[i][:])
            println()
    end
    println()
    println()
    println()

    # Comparing Bottom Row
    println("Comparing Bottom")
    for i in 1:matrixObject.numBlocks
            start_index = (i - 1) * matrixObject.blockSize + 1
            end_index = i * matrixObject.blockSize
            println(matrix[end-matrixObject.blockSize+1:end, start_index:end_index])
            println(bottom[i][:])
            println()
    end
    println()
    println()
    println()
end
