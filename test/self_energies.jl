# this file tests the self-energies module
# TODO fix this unit test
module TestSelfEnergies
using QuantumTransport
using Test
using LinearAlgebra

"""
    testElectrodes()

Test function for the `Electrodes` module. It defines the parameters `p` and `A`, creates an `ElectrodesArray`, generates the self-energies `Σₖs`, and checks if `Σₖs` is empty.

"""
function testElectrodes(dict)

    ElectrodesArray = [
            Electrode([-1,0],[0,p.ny],[0,p.nz],p.ny*p.nz,"-x",p.electrodeMaterial,A);
            Electrode([p.nx,p.nx+1],[0,p.ny],[0,p.nz],p.ny*p.nz,"+x",p.electrodeMaterial,A)
    ]
    dict = Dict()

    # Generate the self-energies Σₖs
    Σₖs = genΣₖs(dict,ElectrodesArray) 

    # Print a message to check if Σₖs is empty
    println("Test is checking whether Σₖs is empty")
        @test !isnothing(Σₖs)
    end

# runparams defined in QuantumTransport
testElectrodes(runparams)

end # module