function makeElectrodeH(p::Dict, ElectrodeInfo::Electrode, edge_NNs::Vector{Hopping})
    kfilter = [0; 1; 1] # to ensure that self-energy in x is Γ centered
    function H(k::Vector{Float64})
        # hamiltonian describing the edges
        #Hₑ = spzeros(ComplexF64, 2*p.nsite*p.norb*ElectrodeInfo.n, 2*p.nsite*p.norb*ElectrodeInfo.n)
        rows = Int[]
        cols = Int[]
        elements = ComplexF64[]
        for NN in edge_NNs
            Δϕ = exp(im * (kfilter .* k[1:3]) ⋅ (p["A"] * NN.N))
            #Δϕ = exp(im*k⋅(p.SLa₁*NN.N[1] + p.SLa₂*NN.N[2] + p.SLa₃*NN.N[3]))
            for i = 1:2
                for j = 1:2
                    push!(rows, 2 * NN.b + i - 2)
                    push!(cols, 2 * NN.a + j - 2)
                    push!(elements, copy(NN.t[i, j] * Δϕ))
                end
            end
            #Hₑ[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]*Δϕ
            #Hₑ[2*NN.b  , 2*NN.a-1] += NN.t[2,1]*Δϕ
            #Hₑ[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]*Δϕ
            #Hₑ[2*NN.b  , 2*NN.a  ] += NN.t[2,2]*Δϕ
        end
        return sparse(rows, cols, elements)
    end
    return H
end

# Defines a cᵦ†cₐ term 
function zeeman(Bvals::Vector{Vector{Float64}}, p::Dict, ElectrodeInfo::Electrode)
    # only defined for S-like orbitals with lz = 0
    N = ElectrodeInfo.n * p["nsite"] * p["norb"] * p["nspin"]
    #N = p.n*p.nsite*p.norb*2
    #for i = 1:N
    zeeman = spzeros(ComplexF64, N, N)
    C = ħ / (2 * m₀) #sans q factor -> eV
    #Bxvals = Bvals[:][1]
    for ax = 1:3
        BiVals = [B[ax] for B in Bvals]
        zeeman .+= 2 * C * Diagonal(BiVals) ⊗ I(p["norb"]) ⊗ σ[ax]
    end
    return sparse(zeeman)
end


function HcontactGen(p::Dict, NNs::Vector{Hopping}, ElectrodeInfo::Electrode)
    CedgeNNs = Hopping[]
    edgeNNs = Hopping[]
    LedgeNNs = Hopping[]
    RedgeNNs = Hopping[]
    nextLayerNNs = Hopping[]
    N = ElectrodeInfo.n * p["nsite"] * p["norb"] * 2
    Rvals = RvalsGen(p, ElectrodeInfo)
    Bfield = ElectrodeInfo.A.(Rvals)

    if (p["electrodeMagnetization"] == true)
        Hᵦ = zeeman(Bfield, p, ElectrodeInfo)
    else
        Hᵦ = 0I(N)
    end
    H₀ = spzeros(ComplexF64, N, N) .+ Hᵦ
    #NNs = genNNs(p,ElectrodeInfo)
    NNs = pruneHoppings(NNs, p["prune"]) # cuts off the relevant matrix elements to make thin film
    for NN in NNs
        if (NN.N ⋅ [1; 0; 0] > 0)
            push!(edgeNNs, deepcopy(NN))
            push!(RedgeNNs, deepcopy(NN))
        elseif (NN.N ⋅ [1; 0; 0] < 0)
            push!(edgeNNs, deepcopy(NN))
            push!(LedgeNNs, deepcopy(NN))
        elseif (NN.edge == true)
            push!(edgeNNs, deepcopy(NN))
            push!(CedgeNNs, deepcopy(NN))
        else
            H₀[2*NN.b-1, 2*NN.a-1] += NN.t[1, 1]
            H₀[2*NN.b, 2*NN.a-1] += NN.t[2, 1]
            H₀[2*NN.b-1, 2*NN.a] += NN.t[1, 2]
            H₀[2*NN.b, 2*NN.a] += NN.t[2, 2]
        end
    end
    Hₑ = makeElectrodeH(p, ElectrodeInfo, edgeNNs)
    Hc = makeElectrodeH(p, ElectrodeInfo, CedgeNNs)
    Hₗ = makeElectrodeH(p, ElectrodeInfo, LedgeNNs)
    Hᵣ = makeElectrodeH(p, ElectrodeInfo, RedgeNNs)
    function Hslab(k::Vector{Float64})
        Hcenter = Hc(k)
        if (size(Hcenter) == (0, 0))
            return H₀
        elseif (size(H₀) == (0, 0))
            return Hcenter
        else
            return Hc(k) .+ H₀
        end
    end
    function H(k::Vector{Float64})
        return Hₑ(k) .+ H₀
    end
    return Hslab, Hₗ, Hᵣ, H
end

function RvalsGen(p::Dict, ei::Electrode)
    nsite = p["nsite"]
    N = ei.n * nsite
    R = Vector{Vector{Float64}}(undef, N)
    ix = -1
    if (ei.connectfrom == "-x")
        ix = -1
    elseif (ei.connectfrom == "+x")
        ix = p["nx"] + 1
    else
        return 0 #something is busted
    end
    # offset to fix the 
    iRoffset = Int(0 + 0 + ix * nsite + ei.yrange[1] * p["nx"] * nsite + ei.zrange[1] * p["ny"] * p["nx"] * nsite)
    nx = 1
    for iy = ei.yrange[1]:(ei.yrange[2]-1)
        for iz = ei.zrange[1]:(ei.zrange[2]-1)
            for isite = 0:(nsite-1)
                iR = Int(1 + isite + nsite * (ix - ix + nx * ((iy - ei.yrange[1]) + (iz - ei.zrange[1]) * p["ny"])))
                #iR = Int(1 - iRoffset + isite + ix*p.nsite + iy*p.nx*p.nsite + iz*p.ny*p.nx*p.nsite)
                #println("$ix $iy $iz $isite iR $iR")
                ivec = Int.([ix, iy, iz, isite])
                Rval = xyztor(p, ivec)
                R[iR] = deepcopy(Rval)
            end
        end
    end
    return R # needs ⊗I(p.norb)⊗I(2) for full (spinful) hilbert space
end

