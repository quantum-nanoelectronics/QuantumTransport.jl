module CommonModule
# Do not perform exports here as this may cause compilation issues and clutter the namespace
# Instead, perform exports in the module that imports this module by doing:
# using .CommonModule: ⊗, τ₁, ...

using LinearAlgebra

include("Functions.jl")
include("Data.jl")

# export all symbols

end