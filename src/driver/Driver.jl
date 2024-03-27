using QuantumTransport
using PyPlot

function genBZ(p::Dict,nx::Int=0, ny::Int=100, nz::Int=100) # only works for cubic lattice
    # nx, ny, and nz specifically refer to # of points in IBZ
    kpoints = Vector{Float64}[]
    kindices = Vector{Int}[]
    kweights = Float64[]
    X1 = p["kdict"]["X₁"];
    X2 = p["kdict"]["X₂"];
    X3 = p["kdict"]["X₃"];
    if(nx != 0)
        kxs = collect(LinRange(-X1[1],X1[1],2*nx+1));
    else
        kxs = [0]
    end
    if(ny != 0)
        kys = collect(LinRange(-X2[2],X2[2],2*ny+1));
    else
        kys = [0]
    end
    if(nz != 0)
        kzs = collect(LinRange(-X3[3],X3[3],2*nz+1));
    else
        kzs = [0]
    end
    function divFixNaN(a::Int,b::Int) # for this particular instance, n/0 represents a Γ-centred sampling @ k = 0. 
            if(b==0)
                    return 0
            else
                    return a/b
            end
    end
    for ix = -nx:nx
        for iy = -ny:ny
            for iz = -nz:nz
                kindex = [iy + ny + 1; iz + nz + 1]
                k = divFixNaN(ix,nx)*X1 + divFixNaN(iy,ny)*X2 + divFixNaN(iz,nz)*X3
                kweight = 1
                if(abs(ix) == nx)
                    kweight *= 1/2
                end
                if(abs(iy) == ny)
                    kweight *= 1/2
                end
                if(abs(iz) == nz)
                    kweight *= 1/2
                end
                push!(kpoints,k)
                push!(kindices,kindex)
                push!(kweights,kweight)
            end
        end
    end
    ksum = sum([w for w in kweights])
    kweights = (1/ksum).*kweights
    return kpoints, kweights, kindices, kxs, kys, kzs
end

function SaveFigure(fig,path,name="",type=".svg")
	fig.savefig(path*"/"*name*type)
	PyPlot.close(fig)
        GC.gc()
end

function pyplotHeatmap(x,y,z,xlab="",ylab="",name="",cmap= :nipy_spectral,save=false, path="")
	fig, ax = PyPlot.subplots();
	#=
        #dx = maximum(x)-minimum(x); dy = maximum(y)-minimum(y)
	C = 500
	width = C
	height = C*AR
	#surf = PyPlot.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap, size=(width,height))
	show(size(z))
	show(size(x))
	show(size(y))
        #surf = Plots.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, c = cmap)
        #surf = Plots.heatmap(x,y,z, xlabel=xlab, ylabel=ylab, title=name, cmap = cmap, 
        #                      norm=matplotlib[:colors][:SymLogNorm](1e-2),
        #                     size=(width,height), margin=5mm, zscale=:log10)
        =#
        min=10^-20
        vmin = maximum([min,minimum(z)])
        surf = ax[:imshow](z.+vmin, cmap=cmap, norm=matplotlib[:colors][:LogNorm](vmin=10^0, vmax=10^-1), extent= [minimum(x), maximum(x), minimum(y), maximum(y)])
        ax[:set_xlabel](xlab)
        ax[:set_ylabel](ylab)
        PyPlot.colorbar(surf, label=name)
        if(save)
            SaveFigure(fig,path,name)
            #fig.savefig(path*name*".svg")
            #close(fig)
        end
        #heatmap!
        #xlabel
	#gui(surf)
        return surf
end

function main(p::Dict, A::Function)
    nested_params = generateParams(p)
    emptyElectrode = Electrode([p["nx"],p["nx"]+1],[0,p["ny"]],[0,p["nz"]],p["ny"]*p["nz"],"+x",p["electrodeMaterial"],A)
    #electroces needs to call genNNs/genH itself (not really sure if there's a way to do this unless we pass in a nested array)
    ElectrodesArray = [
            Electrode([-1,0],[0,p["ny"]],[0,p["nz"]],p["ny"]*p["nz"],"-x",p["electrodeMaterial"],A);
            Electrode([p["nx"],p["nx"]+1],[0,p["ny"]],[0,p["nz"]],p["ny"]*p["nz"],"+x",p["electrodeMaterial"],A)
    ]
    p["prune"] = union(["x"], p["prune"])
    p["verbose"] = false
    p["nelectrodes"] = size(ElectrodesArray)[1]

    NNs = genNNs(nested_params["hoppings"])
    NNs = pruneHoppings(NNs, p["prune"])
    H₀, edge_NNs = nnHoppingMat(NNs, nested_params["matrix"])
    returnvals = []
    H = genH(nested_params["hamiltonian"], H₀, edge_NNs, returnvals)

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
    #parallelk = false
    if(p["nk"] > 0)
        S = 0.15 # scale for k-map
    else
        S = 1
    end
    #println("parallelk = $parallelk, negf_params.prune = $(negf_params.prune)")
    Operators = [I(p["nx"]*p["ny"]*p["nz"]*p["norb"]*2), γ⁵]
    TofE, Tmap = totalT(genScatteredT, kindices, S .* kgrid, kweights, p["E_samples"], p["E_samples"][1], parallelk, Operators)
    TofE = S^2*TofE
    #print(Tmap)
    figh = pyplotHeatmap(S*kys/(π/p["a"]),S*kzs/(π/p["a"]),Tmap',"ky (π/a)","kz (π/a)","T(ky,kz)",:nipy_spectral, p["savedata"], p["path"])
    if("tplot" ∈ p["returnvals"])
        push!(returnvals,figh)
    end
end