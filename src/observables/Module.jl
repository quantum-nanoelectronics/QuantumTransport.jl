# renamed to match other module names
# also, more generic names can conflict with functions/data
module ObservablesModule
export genDOS

using LinearAlgebra
using SparseArrays

include("DOS.jl")

end



