module NEGFModule

export NEGF_prep, totalT, DOS, siteDOS, sitePDOS

include("../hoppings/createHoppings.jl")
include("../common/Module.jl")
using .CommonModule: ⊗, ħ, q, eV

using LinearAlgebra
using SparseArrays
using Distributed
using Random

include("NEGF.jl")

end
