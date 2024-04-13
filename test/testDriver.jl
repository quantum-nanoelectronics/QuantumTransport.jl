module testDriver

using QuantumTransport
using LinearAlgebra
include("../data-input/runs.jl")

using Test

include("InBi.jl")

# TODO Remove include statements from /src folder in test files like this one (except do it for all test files)
include("../src/driver/Driver.jl")


function DriverTest(pNNs::NamedTuple, A::Function)
    dict = Dict()
    
    for (key,value) in zip(keys(pNNs), pNNs)
        dict[string(key)] = value
    end
    # println(dict)

    main(dict, A)

    return true
end

function newDriverTest(params)
    main(params)
    return true
end

# p, p1, p2, p3, A = InBi.generateParams()
# @test DriverTest(p, A)

@test newDriverTest(runparams)


end