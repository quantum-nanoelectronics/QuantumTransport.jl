function unitcell(p::Dict)
    G = transpose(2*π*inv(A))
    function kdictGen(A)
        G = transpose(2*π*inv(A))
        kdict = Dict(
            "Γ" => G*[0; 0; 0],
            "A" => G*[1/2; 1/2; 1/2],
            "M" => G*[1/2; 1/2; 0],
            "Z" => G*[0; 0; 1/2],
            #"-Z" => B*[0; 0; -1/2],
            "Y" => G*[0; 1/2; 0],
            #"-Y" => B*[0; -1/2; 0],
            "X" => G*[1/2; 0; 0],
            #"-X" => B*[-1/2; 0; 0],
        )
        return kdict
    end

    NNs = genNNs(p)
    NNs = pruneHoppings(NNs, p["prune"])
    H₀, edge_NNs = nnHoppingMat(NNs, p)
    H = genH(p, p["A_field"], H₀, edge_NNs)

    #observables = Dict(
    #    :bandstructure => () -> getbands(p, p["klist"], kdictGen(p["A"]), p["numInterpolations"], H),
    #    :DOS => () -> genDOS(p, H, 0) # TODO implement DOS
    #)

    
    if haskey(p, "save")
        if :DOS ∈ p["save"]
            DOS = nothing
            if typeof(H) <: Function
                # H is a function of k, and we will integrate over the brillouin zone
                energies = Float64[]
                kpts, dk = genBZ(G, p["nkvec"])
                if typeof(H_Γ) <: Matrix
                    DOS = genDOS(p::Dict, H, type::Symbol, η::Float64)
                    
                else
                    for k ∈ kpts
                        append!(energies, eigvals(Array(H(k))))
                    end
                end
            else
                @warn "Other Hamiltonian types not implemented for DOS"
            end
            DOS_vals = DOS.(p["E_samples"])
            save_data_formatted("ℝ→ℝ", p["path"], "DOS.csv", ["DOS (1/eV⋅nm³)","E (eV)"], [p["E_samples"],TofE]; flip_axes=true, title="Transmission")
            println("DOS: ", DOS_vals)
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

function genBZ(G::Matrix, nkvec::Vector{Int})
    kpoints = Vector{Float64}[]
    nG1 = nkvec[1]; nG2 = nkvec2[2]; nG3 = nkvec3[3];
    dk = det(G)/(nG1*nG2*nG3)
    for ix ∈ 1:nG1
        for iy ∈ 1:nG2
            for iz in 1:nG3
                k = G*[ix/nG1; iy/nG2; iz/nG3]
                push!(kpoints,k)
            end
        end
    end
    return kpoints, weight
end