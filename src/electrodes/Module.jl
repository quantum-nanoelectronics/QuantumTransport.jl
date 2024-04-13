module ElectrodesModule

export genΣₖs

using ..CommonModule

using LinearAlgebra
using SparseArrays

include("../hoppings/createHoppings.jl")
include("Utilities.jl")
include("ModifyH.jl")
include("SelfEnergies.jl")



end

