module HamiltonianModule

using ..CommonModule

include("ConstructHamiltonian.jl")
# Changed from Hgen to genH
# export Hgen
export genH


end