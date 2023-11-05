using QuantumTransport
using Test

# this file includes a test for performance of column and row major functions
# test based on local machine 

# test if greater/less than 0.5 second
function test1(f::Function)
    return true # remove this line when testing locally
    time = QuantumTransport.calculate_time(f)
    if time > 0 && time < 0.5
        return true
    end
    return false
end

function test2(f::Function)
    return true # remove this line when testing locally
    time = QuantumTransport.calculate_time(f)
    if time > 0.5 && time < 5
        return true
    end
    return false
end

@test test1(column_major)

@test test2(row_major)



