function unitcell(p::Dict, A::Function)
    function kdictGen(A)
        B = transpose(2*π*inv(A))
        kdict = Dict(
            "Γ" => B*[0; 0; 0],
            "A" => B*[1/2; 1/2; 1/2],
            "M" => B*[1/2; 1/2; 0],
            "Z" => B*[0; 0; 1/2],
            "-Z" => B*[0; 0; -1/2],
            "X₂" => B*[0; 1/2; 0],
            "-X₂" => B*[0; -1/2; 0],
            "X₁" => B*[1/2; 0; 0],
            "-X₁" => B*[-1/2; 0; 0],
            "X₃" => B*[0; 0; 1/2]
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
        for f in p["save"]
            if haskey(observables, f)
                observables[f]()
            end
        end
    end
    
    
end