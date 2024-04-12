using LinearAlgebra

include("Transport.jl")

function main(p::Dict)
    # println("p Dictionary: ")
    println(p)
    # println()
    # println("A Function: ")
    # println(A)

    # # println(q, ħ)


    # transport(p, A)


    function A(R::Vector{Float64})
        return [0.0, 0.0, 0.0]
    end
    NEGF_Transport_1D(p, A)



    if haskey(p, "unitcell")
        println("=============== Running unitcell ===============")

    end
    if haskey(p, "transport")
        println("=============== Running transport ===============")
        NEGF_Transport_1D(p["transport"], A)
       
    end
    if haskey(p, "supercell")
        println("=============== Running supercell ===============")
    end




    

end
