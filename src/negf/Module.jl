module NEGFModule

export NEGF_prep, totalT, DOS, siteDOS, sitePDOS

using ..CommonModule

include("../hoppings/createHoppings.jl")

using LinearAlgebra
using SparseArrays
using Distributed
using Random

include("NEGF.jl")

end
