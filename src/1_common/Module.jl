# currently, QuantumTransport will export all symbols from common here. Fix this later.

# Include all functions and data in this module if it needs to be used in another submodule in this package.

module CommonModule # the _ in the directory is needed because this needs to be compiled first to be used within other package modules

using LinearAlgebra

include("Functions.jl")
include("Data.jl")
include("Structs.jl")

export ⊗ # for some reason, this needs its own export

for n in names(@__MODULE__; all=true)
    if Base.isidentifier(n) && n ∉ (:eval, :include, Symbol(@__MODULE__))
        @eval export $n
    end
end

# commented out so that updates to input files propogate immediately for visualization on frontend
# include(joinpath(INPUT_DIR, AllInputs.jl))
# export runparams

end
