# currently, QuantumTransport will export all symbols from common here. Fix this later.

module CommonModule

using LinearAlgebra

include("Functions.jl")
include("Data.jl")
include("Structs.jl")


include("../../data-input/materials.jl")

# for some reason, this needs its own export
export ⊗

#exports all units
for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
                #    println("Exported: ", n)
               end
end

end




