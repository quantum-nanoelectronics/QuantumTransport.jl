module CommonModule
# Do not perform exports here as this may cause compilation issues and clutter the namespace
# Instead, perform exports in the module that imports this module by doing:
# using .CommonModule: ⊗, τ₁, ...

using LinearAlgebra

include("Functions.jl")
include("Data.jl")

export ⊗

#exports all units
for n in names(@__MODULE__; all=true)
               if Base.isidentifier(n) && n ∉ (Symbol(@__MODULE__), :eval, :include)
                   @eval export $n
                   println("Exported: ", n)
               end
end

end


