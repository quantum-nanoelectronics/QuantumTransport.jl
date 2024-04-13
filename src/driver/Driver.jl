using LinearAlgebra


# TODO temp
function A_Function(R::Vector{Float64})
    return [0.0, 0.0, 0.0]
end 

include("Transport.jl")

function main(p::Dict, A::Function = A_Function)
    # println("p Dictionary: ")
    # println(p)
    
    # TODO move this in the if statements
    NEGF_Transport_1D(p, A)

    if haskey(p, "unitcell")
        println("=============== Running unitcell ===============")

    end
    if haskey(p, "transport")
        println("=============== Running transport ===============")
        # NEGF_Transport_1D(p["transport"], A) # TODO fix this
       
    end
    if haskey(p, "supercell")
        println("=============== Running supercell ===============")
    end

end
