module testHoppings # TODO fix this unit test

using QuantumTransport
using LinearAlgebra
using Test

function HoppingsTest(dict) 
    NNs = genNN(dict)
    checkNNs(NNs, dict)
    NNs = pruneHoppings(NNs, []) #p needs a prune value
    H₀, edge_NNs = nnHoppingMat(NNs, dict)
    H = testHGen(dict, A, H₀, edge_NNs)
end

function genNN(pNNs) 
    NNs = genNNs(pNNs)

    @assert length(NNs) == 7
    @assert NNs[2].N[1] == -1 && NNs[2].edge == true #-x
    @assert NNs[3].N[1] == 1 && NNs[3].edge == true #+x
    @assert NNs[4].N[2] == -1 && NNs[4].edge == true #-y
    @assert NNs[5].N[2] == 1 && NNs[5].edge == true #+y
    @assert NNs[6].N[3] == -1 && NNs[6].edge == true #-z
    @assert NNs[7].N[3] == 1 && NNs[7].edge == true #+z
    return NNs
end

function checkNNs(NNs, pNNs) 
    newNNs = genNNs(pNNs)
    @testset verbose = true "Hopping Parameters[$i]" for i in range(1,7)
        @test newNNs[i].a == NNs[i].a 
        @test newNNs[i].b == NNs[i].b
        @test newNNs[i].ia == NNs[i].ia 
        @test newNNs[i].ib == NNs[i].ib
        @test newNNs[i].ra == NNs[i].ra 
        @test newNNs[i].rb == NNs[i].rb
        @test newNNs[i].t == NNs[i].t
        @test newNNs[i].edge == NNs[i].edge
        @test newNNs[i].N == NNs[i].N
    end
end

function testHGen(p, A, H₀, edge_NNs)
    H = genH(p, A, H₀, edge_NNs, [])
    sparseArray = H([0;0;0])
    #up-spin, down-spin
    @assert sparseArray.m == 2 && sparseArray.n == 2
end

⊗(A,B) = kron(A,B)

# runparams defined in QuantumTransport
HoppingsTest(runparams)
end

