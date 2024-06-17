using LinearAlgebra

# this statement may cause issues because it is using the package that this file is being included in
# not sure if including the package here is good practice but the alternative is to include it in the testDriver.jl file like so: include("../src/driver/Driver.jl")
using QuantumTransport

include("Transport.jl")
include("Unitcell.jl")
include("Supercell.jl")

function main(p::Dict)
    if haskey(p, "transport")
        println("=============== Running transport ===============")
        @time NEGF_Transport_1D(p["transport"])
    end
    if haskey(p, "unitcell")
        println("=============== Running unitcell ===============")
        unitcell(p["unitcell"])
    end
    if haskey(p, "supercell")
        println("=============== Running supercell ===============")
        supercell(p["supercell"])
    end

end
