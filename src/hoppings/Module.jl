module HoppingsModule

include("../common/Module.jl")

using .CommonModule: ħ, m₀, eV# , σ

using LinearAlgebra

include("createHoppings.jl")
export genNNs, nnHoppingMat, pruneHoppings

end