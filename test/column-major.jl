"""
    test1(f::Function)

Test function to check if the execution time of the function `f` is greater than 0 and less than 0.5 seconds.

# Arguments
- `f::Function`: The function to be tested.

# Returns
- `true` if the execution time of `f` is greater than 0 and less than 0.5 seconds, `false` otherwise.
"""
function test1(f::Function, A::Matrix{Int64})
    
    time = calculate_time(f, A)
    if time > 0 && time < 0.5
        return true
    end
    return false
end



"""
    test2(f::Function)

Test function to check if the execution time of the function `f` is greater than 0.5 and less than 5 seconds.

# Arguments
- `f::Function`: The function to be tested.

# Returns
- `true` if the execution time of `f` is greater than 0.5 and less than 5 seconds, `false` otherwise.
"""
function test2(f::Function, A::Matrix{Int64})
    time = calculate_time(f, A)
    if time > 0.5 && time < 5
        return true
    end
    return false
end

const sampleMatrixA = reshape(collect(1:100000000),10000,10000)

@test test1(column_major, sampleMatrixA)

@test test2(row_major, sampleMatrixA)
