module TestColumnMajor
using QuantumTransport
using Test
println("\033[1mRunning Column Major Tests\033[0m")

"""
    timeColumnMajor(f::Function, A::Matrix{Int64})

This function takes a function `f` and a matrix `A` as input and returns the time taken to execute the function `f` on matrix `A`.

# Arguments
- `f::Function`: The function to be executed on matrix `A`.
- `A::Matrix{Int64}`: The matrix on which the function `f` is executed.

# Returns
- `Float64`: The time taken to execute the function `f` on matrix `A`.

"""
function timeColumnMajor(f::Function, A::Matrix{Int64})
    return calculate_time(f, A)
end

"""
    timeRowMajor(f::Function, A::Matrix{Int64})

Calculate the time taken to execute the function `f` on the matrix `A` using row-major order.

# Arguments
- `f::Function`: The function to be executed on the matrix `A`.
- `A::Matrix{Int64}`: The matrix on which the function `f` will be executed.

# Returns
- The time taken to execute the function `f` on the matrix `A`.

"""
function timeRowMajor(f::Function, A::Matrix{Int64})
    return calculate_time(f, A)
end

sampleMatrixA = reshape(collect(1:100000000),10000,10000)

@test timeColumnMajor(column_major, sampleMatrixA) < timeRowMajor(row_major, sampleMatrixA)

end
