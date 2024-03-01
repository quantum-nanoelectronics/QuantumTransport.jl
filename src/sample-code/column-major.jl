# Sample code to show how to use tests in Julia.
# this file contains functions that measure runtime of julia's row/column major implementations
# Julia stores 2D arrays in column major format, so the corresponding function will access memory faster
# More cache misses will occur in the row_major function

using BenchmarkTools

function column_major(A::Matrix{Int64})
    s = 0
    for j in 1:size(A, 2)
        for i in 1:size(A, 1)
            s += A[i,j]
        end
    end
    return s
end

function row_major(A::Matrix{Int64})
    s = 0
    for i in 1:size(A, 1)
        for j in 1:size(A, 2)
            s += A[i,j]
        end
    end
    return s
end


function calculate_time(f::Function, A::Matrix{Int64})
    return @belapsed $f($A)
end
