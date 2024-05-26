# this file tests the self-energies module
module TestSelfEnergies
using QuantumTransport
using Test
using LinearAlgebra

# Sample A function
A(R::Vector{Float64}) = [0.0, 0.0, 0.0]

include(joinpath(INPUT_DIR, "AllInputs.jl"))

"""
    testElectrodes(p, A)

This function tests the electrodes for the quantum transport simulation.

# Arguments
- `p`: The parameters for the simulation.
- `A`: The input data for the electrodes.

# Returns
- None

"""
function testElectrodes(p, A)
    Electrodes = [
            Electrode([-1,0],[0,p["ny"]],[0,p["nz"]],p["ny"]*p["nz"],"-x",p["electrodeMaterial"],A);
            Electrode([p["nx"],p["nx"]+1],[0,p["ny"]],[0,p["nz"]],p["ny"]*p["nz"],"+x",p["electrodeMaterial"],A)
    ]

    p["nelectrodes"] = size(Electrodes)[1]

    
    NNs = genNNs(p)
    NNs = pruneHoppings(NNs, p["prune"])
    
    # Generate the self-energies Σₖs
    Σₖs = genΣₖs(p, Electrodes)

    # Print a message to check if Σₖs is empty
    println("Test is checking whether Σₖs is empty")
        @test !isnothing(Σₖs)
    end

# runparams defined in QuantumTransport
testElectrodes(runparams["transport"], A)

end # module