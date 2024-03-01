using LinearAlgebra
using SparseArrays

# include("Structs.jl")
# include("Utilities.jl")

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

# Defines a cᵦ†cₐ term 
function zeeman(Bvals::Vector{Vector{Float64}},  p::NamedTuple, ElectrodeInfo::Electrode)
	# only defined for S-like orbitals with lz = 0
	N = ElectrodeInfo.n*p.nsite*p.norb*2
	#N = p.n*p.nsite*p.norb*2
	#for i = 1:N
	zeeman = spzeros(ComplexF64, N, N)
	C = ħ/(2*m₀) #sans q factor -> eV
	#Bxvals = Bvals[:][1]
	for ax = 1:3
		BiVals = [B[ax] for B in Bvals]
		zeeman .+= 2*C*Diagonal(BiVals)⊗I(p.norb)⊗σ[ax]
	end
	return sparse(zeeman)
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
                return Hc(k).+H₀
            end
	end	
	function H(k::Vector{Float64})
		return Hₑ(k).+H₀	
	end	
	return Hslab, Hₗ, Hᵣ, H
end

# Generate H₀ and make a list of edge bonds for generating H(k)
function nnHoppingMat(NNs,p)
	N = p.n*p.nsite*p.norb
	H = zeros(ComplexF64,2*N,2*N)
	edgeNNs = Any[]
	for NN in NNs
		if(NN.edge == true)
			push!(edgeNNs,deepcopy(NN))
		else
			H[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]
			H[2*NN.b  , 2*NN.a-1] += NN.t[2,1]
			H[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]
			H[2*NN.b  , 2*NN.a  ] += NN.t[2,2]
		end
	end
	return H, edgeNNs
end


# Add a bond to the list of bonds, given some list of bonds, coefficient in spin basis, index of both sites, and param list
#=function pushHopping!(NNs::Vector, t, ia::Vector{Int}, ib::Vector{Int}, p) 
	a = xyztoi(p,ia); b = xyztoi(p,ib);
	ra = xyztor(p,ia); rb = xyztor(p,ib); r = rb - ra;
	# for hopping term
	NN = deepcopy(Hopping(a,b,ia,ib,ra,rb,r,t, false, [0;0;0],""))
	push!(NNs,NN)
end=#



#=function weylHopping(p::NamedTuple,NNs::Vector{Hopping},ia::Vector{Int})
        iorb = ia[5]
        t = 3*nextsite(iorb)*p.t*(I(2))
	pushHopping!(NNs, t, ia, ia, p)
	for ax = 1:3
		for dir = [-1,1]
			# for weyl term in hamiltonian
			di = zeros(5); di[ax] = dir; di[5] = nextsite(iorb); ib = Int.(ia + di)
			Ra = xyztor(p,ia); Rb = xyztor(p,ib); 
			δ = Rb - Ra
			# implement H = +vf*𝐩⋅𝛔 = -vf𝑖ħ ∇ᵣ⋅σ on finite grid
			t = (-im/2)*dir*p.t*σ[ax]
			if(any(isnan,t)==true)
				throw(DomainError(t, "Something broken in hamiltonian definition! Returning NaN"))
				return
			end
			pushHopping!(NNs, t, ia, ib, p)
			# for normal hopping term in hamiltonian
			ib[5] = iorb; 
			
			t = -(1/2)*nextsite(iorb)*p.t*(I(2))
			pushHopping!(NNs, t, ia, ib, p)
		end
	end
end=#
	