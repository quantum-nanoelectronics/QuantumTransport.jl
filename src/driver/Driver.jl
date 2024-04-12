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


    if haskey(p, "unitcell")
        println("=============== Running unitcell ===============")

    end
    if haskey(p, "transport")
        println("=============== Running transport ===============")
        transport(p, A)
       
    end
    if haskey(p, "supercell")
        println("=============== Running supercell ===============")
    end




    

end
