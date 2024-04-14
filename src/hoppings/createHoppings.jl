# moved xyztoi and xyztor to CommonModule because !pushHoppings uses them

# Generate H₀ and make a list of edge bonds for generating H(k)
function nnHoppingMat(NNs, p)
    H = zeros(ComplexF64, p["n"], p["n"])
    edgeNNs = Any[]
    for NN in NNs
        if (NN.edge == true)
            push!(edgeNNs, deepcopy(NN))
        else
            # at site 1->2 would be a_2,1
            # with site ⊗ spin, would be 
            H[2*NN.b-1, 2*NN.a-1] += NN.t[1, 1]
            H[2*NN.b, 2*NN.a-1] += NN.t[2, 1]
            H[2*NN.b-1, 2*NN.a] += NN.t[1, 2]
            H[2*NN.b, 2*NN.a] += NN.t[2, 2]
        end
    end
    return H, edgeNNs
end


function electrodeParams(p::Dict, ElectrodeInfo::Electrode)
    nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
    ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
    nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
    ep = deepcopy(p)
    ep["nx"] = nx
    ep["ny"] = ny
    ep["nz"] = nz
    ep["n"] = nx * ny * nz
	ep["deviceMaterial"] = ElectrodeInfo.type
    return ep
end


function genNNs(p::Dict, Electrodes::Electrode)
	ep = electrodeParams(p, Electrodes)
	ep["nx"] = 1
    ep["prune"] = filter(x -> x ∉ ["x"], p["prune"])
	return genNNs(ep), ep
end

#p is a dict
function genNNs(p::Dict) # all of the terms in the hamiltonian get added here, get back the relevant bonds
    NNs = Hopping[]
    nx = p["nx"]
	ny = p["ny"]
	nz = p["nz"]
    hoppingType! = p["material_hamiltonian"][p["deviceMaterial"]]
    # loop over each unit cell site in the superlattice
    for iy = 0:(ny-1)
        for iz = 0:(nz-1)
			for ix = 0:(nx-1)
				for isite = 0:(p["nsite"]-1)
					for iorb = 0:(p["norb"]-1)
						ia = copy(Int.([ix, iy, iz, isite, iorb]))
						hoppingType!(p, NNs, ia)
					end
				end
			end
        end
    end
    # now fix the designation for the vectors that hop out of the lattice
    # Will connect them around later using bloch's theorem to generate H(k) function
    for NN in NNs
        ib = [NN.ib[1], NN.ib[2], NN.ib[3]]
        # Δ(ib,ib reflected back into 1st lattice)
        pib = ib - [mod(ib[1], nx), mod(ib[2], ny), mod(ib[3], nz)]
        if (pib ⋅ pib != 0) # if vector is distinctly outside of 1st lattice
            NN.N = Int.([round(pib[1] / (nx)), round(pib[2] / ny), round(pib[3] / nz)])
            NN.b = xyztoi(p, NN.ib, NN.N)
            NN.edge = true
        end
    end
    return NNs
end

function pruneHoppings(NNs, type)
    if ("x" ∈ type)
        deleteat!(NNs, findall(NN -> NN.N[1] != 0, NNs))
    end
    if ("y" ∈ type)
        deleteat!(NNs, findall(NN -> NN.N[2] != 0, NNs))
    end
    if ("z" ∈ type)
        deleteat!(NNs, findall(NN -> NN.N[3] != 0, NNs))
    end
    #display(NNs)
    return NNs
end
