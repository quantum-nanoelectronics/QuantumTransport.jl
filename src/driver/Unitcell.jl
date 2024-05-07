function unitcell(p::Dict, A::Function)
    function kdictGen(A)
        B = transpose(2*π*inv(A))
        kdict = Dict(
            "Γ" => B*[0; 0; 0],
            "A" => B*[1/2; 1/2; 1/2],
            "M" => B*[1/2; 1/2; 0],
            "Z" => B*[0; 0; 1/2],
            #"-Z" => B*[0; 0; -1/2],
            "Y" => B*[0; 1/2; 0],
            "-Y" => B*[0; -1/2; 0],
            "X" => B*[1/2; 0; 0],
            "-X" => B*[-1/2; 0; 0],
        )
        return kdict
    end

    NNs = genNNs(p)
    NNs = pruneHoppings(NNs, p["prune"])
    H₀, edge_NNs = nnHoppingMat(NNs, p)
    H = genH(p, A, H₀, edge_NNs)

    observables = Dict(
        :bandstructure => () -> getbands(p, p["klist"], kdictGen(p["A"]), p["numInterpolations"], H),
        :DOS => () -> genDOS(p, H, 0) # TODO implement DOS
    )

    
    if haskey(p, "save")
        if :DOS ∈ p["save"]
            if DOS_method == :Gᴿ
            Gᴿ(E,k) = inv((E+im*p["η"])*I(p["n"]) - H(k))
            DOS = genDOS
            elseif DOS_method == :integrate_BZ
                println("Test")
            end
        end
        if :bandstructure ∈ p["save"]
            println("Test")
        end
        #=for f in p["save"]
            if haskey(observables, f)
                observables[f]()
            end
        end=#
    end
    
    
end