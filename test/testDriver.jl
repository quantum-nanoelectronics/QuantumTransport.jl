module testDriver

using QuantumTransport
using LinearAlgebra
using Test

#include("InBi.jl")

function DriverTest(pNNs::NamedTuple, A::Function)
    dict = Dict()
    
    for (key,value) in zip(keys(pNNs), pNNs)
        dict[string(key)] = value
    end
    # println(dict)

    main(dict, A)

    return true
end

function printDict(runparams, type)
    println("\033[1m$(type):\033[0m")
    for (key, value) in runparams
        println("  \033[1m$(key):\033[0m $(value)")
    end
end

function newDriverTest(params)
    main(params)
    return true
end


printDict(runparams["transport"], "Transport Parameters")
printDict(runparams["unitcell"], "Unitcell Parameters")
# printDict(runparams["supercell"], "Supercell Parameters")

@test newDriverTest(runparams)



# old implementation (can use for debugging)
# p, p1, p2, p3, A = InBi.generateParams()
# @test DriverTest(p, A)


end