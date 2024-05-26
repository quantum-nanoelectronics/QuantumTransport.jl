module testHoppings

using QuantumTransport
using LinearAlgebra
using Test

include(joinpath(INPUT_DIR, "AllInputs.jl"))

A(R::Vector{Float64}) = [0.0, 0.0, 0.0]

"""
    HoppingsTest(p, A)

This function performs a test on the hoppings in the QuantumTransport module.

# Arguments
- `p::Int`: The number of particles.
- `A::Array{Float64,2}`: The hopping matrix.

# Returns
- `result::Bool`: `true` if the test passes, `false` otherwise.
"""
function HoppingsTest(p, A) 
    NNs = genNNs(p)
    @test !isnothing(NNs)

    NNs = pruneHoppings(NNs, p["prune"])
    @test !isnothing(NNs)

    H₀, edge_NNs = nnHoppingMat(NNs, p)
    H = genH(p, A, H₀, edge_NNs)
    sparseArray = H([0;0;0])
    @test !isnothing(sparseArray.m)
    @test !isnothing(sparseArray.n)
end

# runparams defined in QuantumTransport
HoppingsTest(runparams["transport"], A)
end

