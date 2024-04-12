module HoppingsModule

include("../common/Module.jl")
using .CommonModule

using LinearAlgebra

include("createHoppings.jl")
export genNNs, nnHoppingMat, pruneHoppings

end