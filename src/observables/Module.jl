# renamed to match other module names
# also, more generic names can conflict with functions/data
module ObservablesModule
export genDOS, getbands

using LinearAlgebra
using SparseArrays

using ..InputOutputModule

include("DOS.jl")
include("Bands.jl")

end



