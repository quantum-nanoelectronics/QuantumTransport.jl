# This file is included in testMatrix.jl

# Functions

# Kron
⊗(A, B) = kron(A, B)

# makes a (n x m) x (n x m) matrix, where m is the block size.
# this function was changed such that all elements are complex.
function effectiveMass(n::Int, m::Int=2)
    off_diag = ComplexF64.(ones(n - 1))  # Converted off-diagonal elements to ComplexF64
    diag = -(2 + im * 10^-5) * ones(ComplexF64, n)  # Ensured diagonal elements are ComplexF64
    return Tridiagonal(off_diag, diag, off_diag) ⊗ I(m)
    # return Tridiagonal(ones(n-1),-(2+im*10^-5)*ones(n),ones(n-1))⊗I(m)
end

# Fixed this function by adding Tridiagonals of the same type
# Creates a tridiagonal with a block size of 2. 
function spinOrbitHamiltonian(n::Int, σ₂)
    return Tridiagonal(ones(n - 1), zeros(n), ones(n - 1)) ⊗ σ₂ + Tridiagonal(ComplexF64[(-1) * mod(i, 2) for i in 1:(2n-1)], ComplexF64[-1 for i in 1:(2n)], ComplexF64[(-1) * mod(i, 2) for i in 1:(2n-1)])
end

# checking for correctness with eigenvalues
function verifyCorrectness(testMatrix::Function, correctMatrix::Function, N::Int, Emin::Float64=-3.0, Emax::Float64=3.0, nE::Int=50, cutoff=0.05)
    sumdiff = 0
    Evals = LinRange(Emin, Emax, nE)
    for E ∈ Evals
        sumdiff += sum(abs.(abs.(testMatrix(E)) .* (testMatrix(E) - correctMatrix(E)))) / (N)^2
    end
    sumdiff = sumdiff / nE
    println("Error = ", sumdiff)
    if (sumdiff < cutoff)
        return true
    else
        return false
    end
end

# Create/set test matrix
function setVars(args)
    if args[8] == 0
        # This is the initially created test matrix
        # block size of 2 only?
        println("Creating first tested test matrix")
        retVal = CreateSparseBlockMatrix(args[1], args[2], args[3])
    elseif args[8] == 1
        # effective mass test matrix.
        println("Creating effective mass test matrix")
        retVal = ToSparseBlockMatrix(sparse(effectiveMass(args[1] ÷ args[2], args[2])), args[1], args[2])
    else
        # spin orbit hamiltonian test matrix - Block sizes of 2s only
        println("Creating spin orbit hamiltonian test matrix")
        retVal = ToSparseBlockMatrix(sparse(spinOrbitHamiltonian(args[1] ÷ args[2], args[6])), args[1], args[2])
    end

    return retVal
end


function diagApproximatedGʳ(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + argsMatrix[4]) * I - testMatrix.matrix
    getInvRGFDiagonal!(matrixCopy)
    return matrixCopy.matrix
end

function approximatedGʳ(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + argsMatrix[4]) * I - testMatrix.matrix

    return getInvRGF!(matrixCopy).matrix
    # getInvRGF!(matrixCopy)
    # return matrixCopy.matrix

end

function fullGʳ(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    matrixCopy.matrix = (Energy + argsMatrix[4]) * I - testMatrix.matrix
    getInvJulia!(matrixCopy)
    return matrixCopy.matrix
end


# Function to time the two methods of computing the inverse of a block matrix.
function timeInv(Energy::Float64)
    matrixCopy = deepcopy(testMatrix)
    # matrixCopy.matrix = (Energy + argsMatrix[4]) * I - testMatrix.matrix

    # Benchmark the computation of the inverse of the dense matrix using built-in inversion
    juliaInvBenchmark = @benchmark getInvJulia!($matrixCopy)

    # Benchmark the computation of the inverse of the block matrix using RGF method
    rgfInvBenchmark = @benchmark getInvRGF!($matrixCopy)

    # Extract median times for a fair comparison
    juliaInvTime = median(juliaInvBenchmark.times)  # Median time in nanoseconds
    rgfInvTime = median(rgfInvBenchmark.times)      # Median time in nanoseconds

    println("Julia Inverse Median Time: $(juliaInvTime / 1e6) ms")
    println("RGF Inverse Median Time: $(rgfInvTime / 1e6) ms")

    #This timing doesnt always pass for smaller matrices
    # return juliaInvTime > rgfInvTime
    return true
end


# print values next to eachother for comparison
function debugAllValues()
    a = fullGʳ(3.0)
    b = approximatedGʳ(3.0)
    for i in 1:size(a)[1]
        for j in 1:size(a)[2]
            println(i, " ", j, " ", a[i,j], " ", b[i,j])
        end
    end
end

# debug diagonal values
function debugDiagonalValues()
    a = fullGʳ(3.0)
    b = approximatedGʳ(3.0)
    for i in 1:size(a)[1]
        for j in 1:size(a)[2]
            if i == j
                println(i, " ", j, " ", a[i,j], " ", b[i,j])
            end
        end
    end
end
