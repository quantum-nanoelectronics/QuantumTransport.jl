module ElectrodesModule

export Electrode, genΣₖs

include("../common/Module.jl")

# import .CommonModule: ⊗

# fix the imports and using statements, the rest are in Materials

using .CommonModule: ħ, m₀, eV# , σ

using LinearAlgebra
using SparseArrays

include("../hoppings/createHoppings.jl")
include("Utilities.jl")
include("ModifyH.jl")
include("SelfEnergies.jl")



end






#=
module Electrodes
using LinearAlgebra
using Constants
using Operators
using UsefulFunctions
using SparseArrays
#using Materials
using VectorPotential


export Electrode, genΣₖs


include("../common/Module.jl")

# import .CommonModule: ⊗

# fix the imports and using statements, the rest are in Materials

using .CommonModule: ħ, m₀, eV# , σ

using LinearAlgebra
using SparseArrays

include("../hoppings/createHoppings.jl")
include("Utilities.jl")
include("ModifyH.jl")
include("SelfEnergies.jl")

mutable struct Electrode
	xrange::Vector{Int}
	yrange::Vector{Int}
	zrange::Vector{Int}
	n::Int
	connectfrom::String # direction connecting fThreads.@threadsrom (i.e, "-x","+x",etc)
	type::String # material that the electrode is made of
	A::Function # magnitude of the exchange field
end


#=mutable struct Hopping
	a::Int # orbital/site index 1
	b::Int # orbital/site index 2 with PBC
	ia # index vector of site A
	ib # index vector of site B without PBC
	ra # location of atom a
	rb # location of atom b
	r # radius from a to b
	t  # hopping parameter affiliated with c†₂c₁ in spin basis. (i.e., t*I(2) or t*σ₊ might make sense)
	edge::Bool # does this hop off the edge of the superlattice?
	N # vector describing the [n₁;n₂;n₃]⋅[a₁;a₂;a₃] superlattice unit cell of site ib
	desc::String
end=#

function xyztoi(p,ivec, N::Vector{Int} = [0;0;0]) 
	# indexing 0 to N-1
	# in case the A₂ lattice vector != C*[0;1;0]
	diy = Int(round(p.SLa₂[1]/p.a₁[1]))*N[2]
	#ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4]; iorb = ivec[5]
	ix = mod(ivec[1],p.nx); iy = mod(ivec[2] + diy,p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	#ix = mod(ivec[1],p.nx); iy = mod(ivec[2],p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	return iorb + p.norb*isite + p.nsite*p.norb*ix + p.nsite*p.norb*p.nx*iy + p.nsite*p.norb*p.nx*p.ny*iz + 1
end

# include("Materials.jl")
# Same as above, except returns the corresponding atomic position of each index vector 
# useful for calculating ∫A⋅δR peierls phase
function xyztor(p,ivec)
    ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4];
    δdict = Dict(0 => p["A"]*[0.0; 0.0; 0.0], #In 1 
             1 => p["A"]*[0.5; 0.5; 0.5]) #In 2
             #2 => p["A"]*[0.0; 0.5; 0.8975-0.5], #Bi 1
             #3 => p["A"]*[0.5; 0.0; 1.10248-0.5]) #Bi 2
    δ = δdict[isite]
    R = p["a₁"]*ix + p["a₂"]*iy + p["a₃"]*iz + δ
        #println("ivec = $ivec, Rpos = $R")
        return R
end
function RvalsGen(p::Dict{String,Any},ei::Electrode)
    N = ei.n*p["nsite"]
    R = Vector{Vector{Float64}}(undef,N)
    ix = -1
    if(ei.connectfrom=="-x")
        ix = -1
    elseif(ei.connectfrom=="+x")
        ix = p["nx"]+1
    else
        return 0 #something is busted
    end
    # offset to fix the 
    iRoffset = Int(0 + 0 + ix*p["nsite"] + ei.yrange[1]*p["nx"]*p["nsite"] + ei.zrange[1]*p["ny"]*p["nx"]*p["nsite"])
        nx = 1
        for iy = ei.yrange[1]:(ei.yrange[2]-1)
            for iz = ei.zrange[1]:(ei.zrange[2]-1)
                for isite = 0:(p["nsite"]-1)
                    iR = Int(1 + isite + p["nsite"]*(ix-ix + nx*((iy-ei.yrange[1]) + (iz-ei.zrange[1])*p["ny"])))
                    #iR = Int(1 - iRoffset + isite + ix*p["nsite"] + iy*p["nx"]*p["nsite"] + iz*p["ny"]*p["nx"]*p["nsite"])
                    #println("$ix $iy $iz $isite iR $iR")
                    ivec = Int.([ix,iy,iz,isite])
                    Rval = xyztor(p,ivec)
                    R[iR] = deepcopy(Rval)
                end
            end
    end
    return R # needs ⊗I(p["norb"])⊗I(2) for full (spinful) hilbert space
end
# Defines a cᵦ†cₐ term 
function zeeman(Bvals::Vector{Vector{Float64}},  p::Dict{String,Any}, ElectrodeInfo::Electrode)
    # only defined for S-like orbitals with lz = 0
    N = ElectrodeInfo.n*p["nsite"]*p["norb"]*2
    zeeman = spzeros(ComplexF64, N, N)
    C = p["ħ"]/(2*p["m₀"]) #sans q factor -> eV
    for ax = 1:3
        BiVals = [B[ax] for B in Bvals]
        zeeman .+= 2*C*Diagonal(BiVals)⊗I(p["norb"])⊗σ[ax]
    end
    return sparse(zeeman)
end

#=function genΣₖs(p::NamedTuple,ElectrodeInfo::Vector{Electrode})   
	nE = size(ElectrodeInfo)[1]
	Σks = Vector{Function}(undef,nE)
	for i = 1:nE # define a self-energy function for every electrode attached to device
            if(p.electrodeMaterial=="mtjweyl")
                # for the normal part of the weyl electrodes
                NNs = genNNs(p,ElectrodeInfo[i])
		P = changeBasis(p,ElectrodeInfo[i])
		Hslab, Hₗ, Hᵣ = HcontactGen(p,NNs,ElectrodeInfo[i]) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
                # doing some jank to get the coupling matrices right
                insElectrode = deepcopy(ElectrodeInfo[i])
                insElectrode.type = "wins"
                insNNs = genNNs(p,insElectrode)
		insHslab, insHₗ, insHᵣ = HcontactGen(p,NNs,insElectrode) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
                connectDict = Dict{String,Function}("-x"=>insHᵣ,"+x"=>insHₗ)
                Hᵥ = connectDict[ElectrodeInfo[i].connectfrom]
                Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
		Σks[i] = Σk
            else
                NNs = genNNs(p,ElectrodeInfo[i])
		P = changeBasis(p,ElectrodeInfo[i])
		Hslab, Hₗ, Hᵣ = HcontactGen(p,NNs,ElectrodeInfo[i]) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
                #connectDict = Dict{String,Function}("-x"=>Hᵣ,"+x"=>Hₗ)
                if(ElectrodeInfo[i].connectfrom=="-x")
                    Hᵥ = Hᵣ
                    Σk = k -> Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
                    #Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
                    Σks[i] = deepcopy(Σk)
                    #Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
		    #Σks[i] = Σk
                else
                    Hᵥ = Hₗ
                    Σk = k -> Σgen(p,Hslabk),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
                    #Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
                    Σks[i] = deepcopy(Σk)
                end
                #Hᵥ = connectDict[ElectrodeInfo[i].connectfrom]
            end
            #Σk(k) = Σgen(p,Hslab(k),Hₗ(k).+Hᵣ(k),Hᵥ(k),ElectrodeInfo[i],P)
            #Σk(k) = Σgen(p,Hslab(k),Hₗ(k),Hᵣ(k),ElectrodeInfo[i],P)
	end
	return Σks
end=#

# return a vector of Σ(k) functions which return Σₖ(E) which return a sparse nsite × nsite matrix at a given energy
function genΣₖs(p::NamedTuple,ElectrodeInfo::Vector{Electrode})   
	nE = size(ElectrodeInfo)[1]
	Σks = Vector{Function}(undef,nE)
	for i = 1:nE # define a self-energy function for every electrode attached to device
            NNs = genNNs(p,ElectrodeInfo[i])
            P = changeBasis(p,ElectrodeInfo[i])
            #(Hslab, Hₗ, Hᵣ) = HcontactGen(p,NNs,ElectrodeInfo[i]) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
            Hs = HcontactGen(p,NNs,ElectrodeInfo[i]) 
            # Thus, Hs = (Hslab(k), Hₗ(k), Hᵣ(k))
            # ∃ one contact on left, (nE-1) Contacts on right
            iCinD = Int.(sign(-i+1.5)/2+2.5) # gets the right index for the coupling hamiltonian  
            iCfromD = Int.(sign(i-1.5)/2+2.5) # gets the right index for the coupling hamiltonian  
            if(p.electrodeMaterial=="mtjweyl")
                # doing some jank to get the coupling matrices right
                insElectrode = deepcopy(ElectrodeInfo[i])
                insElectrode.type = "wins"
                insNNs = genNNs(p,insElectrode)
                #println("InsNNs: \n")
                #display(insNNs)
                #insHslab, insHₗ, insHᵣ = HcontactGen(p,NNs,insElectrode) # returns in-electrode matrix, then H_inelectrode(k), H_contact(k) for transverse k
                HsWMTJ = HcontactGen(p,insNNs,insElectrode)
                V  = HsWMTJ[iCinD]
                βₐ = Hs[iCfromD] # away from device
                βₜ = Hs[iCinD] # towards device
            else
                V  = Hs[iCinD]
                βₐ = Hs[iCfromD] # away from device
                βₜ = Hs[iCinD] # towards device
            end
			# println("βₐ11: ", βₐ([0.0]))

			# if (i == 1)
			# 	println("βₐ1: ", size(βₐ))
			# end
            #Σk(k) = Σgen(p,Hs[1](k),Hᵥₑ(k),Hᵥₘ(k),ElectrodeInfo[i],P)

			
            kxes, kweights, kindices = genTetBZ(electrodeParams(p,ElectrodeInfo[1]),1000,0,0)

			println("βₐ11: ", βₐ([0.0]))

            Σk(k) = TΣgen(p,Hs[1](k),βₐ(k),βₜ(k),V(k),ElectrodeInfo[i],P)
            #Σk(k) = Σgen(p,Hs[1](k),βₐ(k),βₜ(k),V(k),ElectrodeInfo[i],P)
            #Σk(k) = ∫Σgen(p,Hs[4], Hs[1](k), βₐ(k), βₜ(k), V(k),ElectrodeInfo[i],P,k,kxes,kweights)
            #Σk(k) = Σgen(p,Hs[1](k),Hs[2](k).+Hs[3](k),Hᵥₘ(k),ElectrodeInfo[i],P)
            Σks[i] = Σk
	end

	return Σks
end

function xyzElectrodeSiteToI(ElectrodeInfo::Electrode,ivec::Vector{Int})
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	   ix = 0*ivec[1]; iy = ivec[2]; iz = ivec[3]
	return  ix + nx*(iy + ny*(iz)) + 1
end

function xyzElectrodeToI(p::NamedTuple, ElectrodeInfo::Electrode,ivec::Vector{Int})
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; iorb = ivec[4]
	return  iorb + p.norb*(ix + nx*(iy + ny*(iz))) + 1
end

function electrodeSiteToDeviceIndex(p::NamedTuple, ElectrodeInfo::Electrode,ivecContact::Vector{Int})
	# now construct an ivec for the site in the device
	if(ElectrodeInfo.connectfrom=="-x")
		ix = 0 # the maximal site in x, edge of electrode
	else # just uhh, presuming that we will only connect in +- x. Can be changed...
		ix = (p.nx-1)
	end
	iy = ivecContact[2]+ElectrodeInfo.yrange[1]
	iz = ivecContact[3]+ElectrodeInfo.zrange[1]
	return ix + p.nx*(iy+p.ny*iz) + 1
end

function changeBasis(p::NamedTuple,ElectrodeInfo::Electrode)
	#nE = ElectrodeInfo.n*p.nsite
	#nD = p.n*p.nsite
	nE = ElectrodeInfo.n
	nD = p.n
        #println("nE = $nE; nD = $nD")
	Psite = spzeros(nD,nE)
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	# only consider sites that will actually touch the slab
	ix = 0
	if(ElectrodeInfo.connectfrom=="-x")
		ix = p.nx-1 # the maximal site in x, edge of electrode
	else # just uhh, presuming that we will only connect in +- x. Can be changed...
		ix = 0
	end
	for iy = 0:(ny-1)
		for iz = 0:(nz-1)
			ivec = [ix,iy,iz,0]
			contactSiteIndex = xyzElectrodeSiteToI(ElectrodeInfo,ivec)
			deviceSiteIndex = electrodeSiteToDeviceIndex(p,ElectrodeInfo,ivec)
			#println("Device site: $deviceSiteIndex, Contact site: $contactSiteIndex")
			Psite[deviceSiteIndex,contactSiteIndex] = 1
		end
	end
        return dropzeros(Psite⊗I(p.nsite*p.norb*2))
        #return sparse(Psite⊗ones(p.nsite,p.nsite)⊗I(p.norb*2))
end
			


function makeElectrodeH(p::NamedTuple,ElectrodeInfo::Electrode,edge_NNs::Vector{Hopping})
	kfilter = [0;1;1] # to ensure that self-energy in x is Γ centered
	function H(k::Vector{Float64})
		# hamiltonian describing the edges
		#Hₑ = spzeros(ComplexF64, 2*p.nsite*p.norb*ElectrodeInfo.n, 2*p.nsite*p.norb*ElectrodeInfo.n)
                rows = Int[]; cols = Int[]; elements = ComplexF64[];
		for NN in edge_NNs
			Δϕ = exp(im*(kfilter.*k)⋅(p.A*NN.N))
			#Δϕ = exp(im*k⋅(p.SLa₁*NN.N[1] + p.SLa₂*NN.N[2] + p.SLa₃*NN.N[3]))
                        for i = 1:2
                            for j = 1:2
                                push!(rows,2*NN.b+i-2);
                                push!(cols,2*NN.a+j-2);
                                push!(elements,copy(NN.t[i,j]*Δϕ))
                            end
                        end
			#Hₑ[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]*Δϕ
			#Hₑ[2*NN.b  , 2*NN.a-1] += NN.t[2,1]*Δϕ
			#Hₑ[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]*Δϕ
			#Hₑ[2*NN.b  , 2*NN.a  ] += NN.t[2,2]*Δϕ
		end
                return sparse(rows,cols,elements)
	end
	return H
end

function electrodeParams(p::NamedTuple,ElectrodeInfo::Electrode)
	nx = Int(abs(ElectrodeInfo.xrange[2] - ElectrodeInfo.xrange[1]))
	ny = Int(abs(ElectrodeInfo.yrange[2] - ElectrodeInfo.yrange[1]))
	nz = Int(abs(ElectrodeInfo.zrange[2] - ElectrodeInfo.zrange[1]))
	ep = (nx = nx, ny = ny, nz = nz, n = nz*nx*ny)
	return merge(p,ep)
end

function HcontactGen(p::NamedTuple,NNs::Vector{Hopping},ElectrodeInfo::Electrode)
	CedgeNNs = Hopping[]
	edgeNNs = Hopping[]
	LedgeNNs = Hopping[]
	RedgeNNs = Hopping[]
	nextLayerNNs = Hopping[]
	N = ElectrodeInfo.n*p.nsite*p.norb*2
	Rvals = RvalsGen(p,ElectrodeInfo)
	Bfield = ElectrodeInfo.A.(Rvals)
	#Bfield = ElectrodeInfo.A.(Rvals)
	if(p.electrodeMagnetization==true)
		Hᵦ = zeeman(Bfield,p,ElectrodeInfo)
	else
		Hᵦ = 0I(N)
	end
	H₀ = spzeros(ComplexF64,N, N) .+ Hᵦ
	#NNs = genNNs(p,ElectrodeInfo)
	NNs = pruneHoppings(NNs,p.prune) # cuts off the relevant matrix elements to make thin film
	for NN in NNs
		if(NN.N⋅[1;0;0] > 0)
			push!(edgeNNs,deepcopy(NN))
			push!(RedgeNNs,deepcopy(NN))
		elseif(NN.N⋅[1;0;0] < 0)
			push!(edgeNNs,deepcopy(NN))
			push!(LedgeNNs,deepcopy(NN))
		elseif(NN.edge==true)
			push!(edgeNNs,deepcopy(NN))
			push!(CedgeNNs,deepcopy(NN))
		else
			H₀[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]
			H₀[2*NN.b  , 2*NN.a-1] += NN.t[2,1]
			H₀[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]
			H₀[2*NN.b  , 2*NN.a  ] += NN.t[2,2]
		end
	end
	Hₑ = makeElectrodeH(p,ElectrodeInfo,edgeNNs)
	Hc = makeElectrodeH(p,ElectrodeInfo,CedgeNNs)
	Hₗ = makeElectrodeH(p,ElectrodeInfo,LedgeNNs)
	Hᵣ = makeElectrodeH(p,ElectrodeInfo,RedgeNNs)
	function Hslab(k::Vector{Float64})
            Hcenter = Hc(k)
            if(size(Hcenter) == (0,0))
                return H₀
            elseif(size(H₀) == (0,0))
                return Hcenter
            else
				println("SIZES", size(Hcenter), size(H₀))
                return Hc(k).+H₀
            end
	end	
	function H(k::Vector{Float64})
		return Hₑ(k).+H₀	
	end	
	return Hslab, Hₗ, Hᵣ, H
end

function ∫Σgen(p::NamedTuple, H::Function, Hslab::SparseMatrixCSC, βₐ::SparseMatrixCSC, βₜ::SparseMatrixCSC, V::SparseMatrixCSC, ElectrodeInfo::Electrode, P::SparseMatrixCSC, k::Vector{Float64}, kxes::Vector{Vector{Float64}}, kweights::Vector{Float64}, cutoff::Float64=10^-5*eV)
    n = ElectrodeInfo.n*p.nsite*p.norb*2
    kvals = [[kₓ[1],k[2],k[3]] for kₓ in kxes]
    function Σ(E::Float64)
        Gₑ = spzeros(ComplexF64,n,n)
        Gₑs = map(k->grInv((E+im*p.η)*I(n) .- H(k)),kvals)
        for i in eachindex(kweights)
                Gₑ .+= kweights[i]*Gₑs[i]
        end
        Gsurf = grInv((E+im*p.η)*I(n) .- Hslab .- βₜ*Gₑ*βₐ)
        Σ_surf = V*Gₑ*V'
        return P*Σ_surf*P'
    end
    return Σ
end

function TΣgen(p::NamedTuple,H::SparseMatrixCSC,βₐ::SparseMatrixCSC, βₜ::SparseMatrixCSC, V::SparseMatrixCSC, ElectrodeInfo::Electrode, P, cutoff::Float64=10^-12*eV)
#function Σgen(p::NamedTuple,H::SparseMatrixCSC,H_coupling::SparseMatrixCSC, Hᵥ::SparseMatrixCSC, ElectrodeInfo::Electrode, P, cutoff::Float64=10^-7*eV)
    n = ElectrodeInfo.n*p.nsite*p.norb*2
    # so H needs to be instantiated and called outside of the loop
    #H_coupling = Hₗ .+ Hᵣ# couples a layer to the infinite on both sides
    #BLAS.set_num_threads(1) # disable linalg multithreading and parallelize over k instead
    function Σ(E::Float64)
		# println("βₐ: ", size(βₐ))
        #Gₑ = grInv((E+im*p.η)*I(n) .- H) # guess for the green's functions in the electrodes
        #Σ_guess = H_coupling*grInv((E+im*p.η)*I(n) .- H .- 0.1*I(n))*Hᵥₘ'
        #Σ_guess = H_coupling*grInv((E+im*p.η)*I(n) .- H .- 0.1*I(n))*H_coupling'
        # converge the self energy 
        error = 1
        # using transfer matrix method described in 
        ω = (E + im*p.η)*I(n)
        t₀ = inv(Array(ω.-H))*βₐ
        #t̃₀ = grInv(ω.-H)*βₜ
        t̃₀ = inv(Array(ω.-H))*βₜ
        #t̃₀ = grInv(ω.-H)*βₜ
        #t₁ = grInv(I(n) .- t₀*t̃₀ .- t̃₀*t₀)*t₀^2
        t̃₋ = I(n); t̃ᵢ = Array(t̃₀)
        t₋ = I(n); tᵢ = Array(t₀)
        Tᵢ = zeros(ComplexF64,n,n)
        Tᵢ .+= t₀
		Π = I(n)
		while error > cutoff
			#println("error = $error, size(T) = $(size(Tᵢ))")
			T₋ = deepcopy(Tᵢ)
			t̃₋ = t̃ᵢ
			t₋ = tᵢ
			#T₋ = deepcopy(Tᵢ)
			#t̃₋ = deepcopy(t̃ᵢ)
			#t₋ = deepcopy(tᵢ)
			Π = Π*t̃₋
			#Π = deepcopy(Π*t̃₋)
			#display(Π)
			#println("")
        	t̃ᵢ = inv(I(n) .- t₋*t̃₋ .- t̃₋*t₋)*t̃₋^2
        	tᵢ = inv(I(n) .- t₋*t̃₋ .- t̃₋*t₋)*t₋^2
                Tᵢ .+= Π*tᵢ
			error = norm(Tᵢ .- T₋)/norm(Tᵢ)	
		end
                #println("Converged, error = $error")
                effH = Array(ω .- H .- βₜ*Tᵢ)
		Gsurf = inv(effH)

		#Gsurf = grInv((E+im*p.η)*I(n) .- H .- conj(V)*Gₑ*βₐ) # guess for the green's functions in the electrodes
        #Gsurf = grInv(
        Σ_surf = V*Gsurf*V'
        #Σ_surf = V*Gsurf*V'
        #Σ_surf = (βₐ)*Gₑ*(V)'
        #Σ_surf = spzeros(ComplexF64,n,n)
        return sparse(P*Σ_surf*P')
    end
    return Σ
end




function genNNs(p,ElectrodeInfo::Electrode) # all of the terms in the hamiltonian get added here, get back the relevant bonds
	n = p.n
	NNs = Hopping[]
	ep = electrodeParams(p,ElectrodeInfo) # electrodeparams
	ix = 0
	nx = 1
	hopping! = hoppingDict[ElectrodeInfo.type]
	for iy = 0:(ep.ny-1)
		for iz = 0:(ep.nz-1)
		    for isite = 0:(ep.nsite-1)
				for iorb = 0:(ep.norb-1)
					ia = (copy(Int.([ix,iy,iz,isite,iorb])));
					hopping!(ep,NNs,ia)
					#if(ElectrodeInfo.type=="weyl")
					#        weylHopping(ep,NNs,ia)
					#end
				end
			end
		end
	end
	# now fix the designation for the vectors that hop out of the lattice
	# Will connect them around later using bloch's theorem to generate H(k) function
	for NN in NNs
		#println("pre hop ($(NN.ia) to $(NN.ib)) = $(NN.a) to $(NN.b)")
		ib = [NN.ib[1],NN.ib[2],NN.ib[3]]
		# Δ(ib,ib reflected back into 1st lattice)
		pib = ib - [mod(ib[1],nx),mod(ib[2],ep.ny),mod(ib[3],ep.nz)]
		#pib = ib - [mod(ib[1],p.nx),mod(ib[2],p.ny),mod(ib[3],p.nz)]
		if(pib⋅pib != 0) # if vector is distinctly outside of 1st lattice
			NN.N = Int.([round(pib[1]/(nx)),round(pib[2]/ep.ny),round(pib[3]/ep.nz)])
			NN.b = xyztoi(ep,NN.ib, NN.N)
			NN.edge = true
			#println("$(NN.N)")
		end
	end
	return NNs
end


function pruneHoppings(NNs, type)
	if("x" ∈ type)
		deleteat!(NNs, findall(NN->NN.N[1]!=0,NNs))
	end
	if("y" ∈ type)
		deleteat!(NNs, findall(NN->NN.N[2]!=0,NNs))
	end
	if("z" ∈ type)
		deleteat!(NNs, findall(NN->NN.N[3]!=0,NNs))
	end
        #display(NNs)
        return NNs
end



end



=#