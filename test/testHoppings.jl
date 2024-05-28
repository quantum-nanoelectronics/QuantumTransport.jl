module testHoppings

using QuantumTransport
using LinearAlgebra
using Test

include(joinpath(INPUT_DIR, "AllInputs.jl"))

"""
    HoppingsTest(p)

This function performs a test on the hoppings in the QuantumTransport module.

# Arguments
- `p::Int`: The number of particles.

# Returns
- `result::Bool`: `true` if the test passes, `false` otherwise.
"""
function HoppingsTest(p) 
    NNs = genNNs(p)
    @test !isnothing(NNs)

    NNs = pruneHoppings(NNs, p["prune"])
    @test !isnothing(NNs)

    H₀, edge_NNs = nnHoppingMat(NNs, p)
    H = genH(p, p["A_field"], H₀, edge_NNs)
    sparseArray = H([0;0;0])
    @test !isnothing(sparseArray)
end

# runparams defined in QuantumTransport
HoppingsTest(runparams["transport"])
end

