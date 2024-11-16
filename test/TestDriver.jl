module testDriver

using QuantumTransport
using LinearAlgebra
using Test

# basically, we only want the "runparams" variable from AllInputs.jl.
include(joinpath(INPUT_DIR, "AllInputs.jl"))

function printDict(runparams, type)
    println("\033[1m$(type):\033[0m")
    for (key, value) in runparams
        println("  \033[1m$(key):\033[0m $(value)")
    end
end

function driverTest(params)
    main(params)
    return true
end

if Base.active_project() != joinpath(BASE_DIR, "Project.toml")
    error("Error: The currently active Julia environment ($(Base.active_project())) is incorrect. Please activate the QuantumTransport environment.")
else
    println("Using the correct environment: $(Base.active_project())")
end

printDict(runparams["transport"], "Transport Parameters")
printDict(runparams["unitcell"], "Unitcell Parameters")
printDict(runparams["supercell"], "Supercell Parameters")

@test driverTest(runparams)


end