module HoppingsModule

using ..CommonModule
using LinearAlgebra

include("createHoppings.jl")
export genNNs, nnHoppingMat, pruneHoppings

end