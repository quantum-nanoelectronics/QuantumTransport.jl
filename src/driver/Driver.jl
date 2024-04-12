using LinearAlgebra

include("Transport.jl")

function main(p::Dict, A::Function)
    println("p Dictionary: ")
    println(p)
    println()
    println("A Function: ")
    println(A)

    println(q, ħ)


    transport(p, A)


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
