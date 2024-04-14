module testDriver

using QuantumTransport
using LinearAlgebra
using Test

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

println("BASE_DIR: ", BASE_DIR)
println("INPUT_DIR: ", INPUT_DIR)
println("OUTPUT_DIR: ", OUTPUT_DIR)


# this will bring runparams into scope
include(joinpath(INPUT_DIR, "AllInputs.jl"))


printDict(runparams["transport"], "Transport Parameters")
# printDict(runparams["unitcell"], "Unitcell Parameters")
# printDict(runparams["supercell"], "Supercell Parameters")

@test newDriverTest(runparams)

end