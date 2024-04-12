function transport(p::Dict, A::Function)
    nested_params = generateParams(p)
    emptyElectrode = Electrode([p["nx"],p["nx"]+1],[0,p["ny"]],[0,p["nz"]],p["ny"]*p["nz"],"+x",p["electrodeMaterial"],A)
    #electroces needs to call genNNs/genH itself (not really sure if there's a way to do this unless we pass in a nested array)
    ElectrodesArray = [
            Electrode([-1,0],[0,p["ny"]],[0,p["nz"]],p["ny"]*p["nz"],"-x",p["electrodeMaterial"],A);
            Electrode([p["nx"],p["nx"]+1],[0,p["ny"]],[0,p["nz"]],p["ny"]*p["nz"],"+x",p["electrodeMaterial"],A)
    ]
    p["prune"] = union(["x"], p["prune"])
    println("prune: ", p["prune"])
    p["verbose"] = false
    p["nelectrodes"] = size(ElectrodesArray)[1]

    NNs = genNNs(nested_params["hoppings"])
    NNs = pruneHoppings(NNs, p["prune"])
    H₀, edge_NNs = nnHoppingMat(NNs, nested_params["matrix"])
    returnvals = []
    H = genH(nested_params["hamiltonian"], A, H₀, edge_NNs, returnvals)

    Σₖs = genΣₖs(p, ElectrodesArray)
    genGᴿ, genT, genA, genScatteredT = NEGF_prep(p, H, Σₖs)

    γ⁵ = I(p["nx"]*p["ny"]*p["nz"])⊗τ₁⊗σ₀
    if(p["mixedDOS"]==true)
        mdE = p["E_samples"][1]*eV; η = 10^-(3.0)
        function plottingGʳ(k::Vector{Float64})
            function Gʳ(E::Float64)
                Σₗ = Σks[1]; Σᵣ = Σks[2]
                return pGrInv((E+im*η)*I(p["nx"]*p["ny"]*p["nz"]*p["norb"]*2) .- H(k) .- Σₗ(k)(E) .- Σᵣ(k)(E),4,false) 
            end
            return Gʳ
        end
        # left and right-handed states
        Operators = [(1/2)*(I(p["nx"]*p["ny"]*p["nz"]*p["norb"]*2).+d*γ⁵)  for d = [-1,+1]]
        DOS = sitePDOS(plot_params,plottingGʳ,Operators, mdE)
        #G = genGʳ([0.1;0.1;0.1])(0.02*eV)
        #Test = DOS([0.1;0.1;0.1]);
        #testDOS(k) = ones(pnx)*(k⋅k)/nm^2;
        nkDOS = 250
        fsfig = mixedDOS(plot_params,DOS,nkDOS,nkDOS)
        if("mixedDOS" ∈ p["returnvals"])
            push!(returnvals,fsfig)
        end
    end

    nkx = p["nk"]*!("x" ∈ p["prune"]); nky = p["nk"]*!("y" ∈ p["prune"]); nkz = p["nk"]*!("z" ∈ p["prune"]);

    kgrid, kweights, kindices, kxs, kys, kzs = genBZ(p,nkx,nky,nkz) # generate surface BZ points
    println("Sweeping transmission over kgrid: $(nkx*2+1), $(nky*2+1), $(nkz*2+1) ...")
    #TofE, Tmap, TmapList = totalT(genT, kindices, 0.3 .* kgrid, kweights, pE_samples, minimum(pE_samples))
    #TofE, Tmap = totalT(genT, kindices, 0.3 .* kgrid, kweights, pE_samples, minimum(pE_samples))

    parallelk = ((nkx+1)*(nky+1)*(nkz+1) > 8)
    # setting parallelk to false for local testing


    if(p["nk"] > 0)
        S = 0.15 # scale for k-map
    else
        S = 1
    end
    #println("parallelk = $parallelk, negf_params.prune = $(negf_params.prune)")
    Operators = [I(p["nx"]*p["ny"]*p["nz"]*p["norb"]*2), γ⁵]

    # testing (remove)
    parallelk = false


    TofE, Tmap = totalT(genScatteredT, kindices, S .* kgrid, kweights, p["E_samples"], p["E_samples"][1], parallelk, Operators)
    TofE = S^2*TofE

    println("TofE: ", TofE)
    #print(Tmap)
    #figh = pyplotHeatmap(S*kys/(π/p["a"]),S*kzs/(π/p["a"]),Tmap',"ky (π/a)","kz (π/a)","T(ky,kz)",:nipy_spectral, p["savedata"], p["path"])
    if("tplot" ∈ p["returnvals"])
        push!(returnvals,figh)
    end

end