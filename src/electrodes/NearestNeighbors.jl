using LinearAlgebra
using SparseArrays

include("Structs.jl")

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


# From GenNNs.jl
# Requires a, A, B, Atom, Atoms
function NNatoms()
	nn1_dist = a/2;
	#AtomNNs = AtomNN[]
	AtomNNs = Array{AtomNN, 1}(undef, (6+8)*4) # each atom will have 6 nn₂s, 8 nn₁s
	n = 1
	for atom = 1:4
		for ax = 1:3
			for dir = [-1,1]
				A1 = Atoms[atom]
				A2 = Atoms[mod((atom+2-1),4)+1]
				δ = zeros(3)
				δ[ax] = dir*(a/2)
				#N = zeros(3)
				
				# lot going on here. This gives distance from imagined periodic
				# atom to the base atom in the unit cell. Take this distance and 
				# multiply by [a1,a2,a3]^-1 (so B) but include (2*π) because
				# this isn't the reciprocal lattice
			
				N = Int.(round.((2*π)^(-1)*B*(δ + A1.R - A2.R)))
				AtomNNs[n] = AtomNN(A1, A2, δ, N)
				n += 1
			end
		end
		for d1 = [-1,1]
			for d2 = [-1,1]
				for d3 = [-1,1]
					dn = d1*d2*d3
					δ = (a/4)*[d1;d2;d3]
					A1 = Atoms[atom]
					A2 = Atoms[mod(atom+dn-1,4)+1]
					N = Int.(round.((2*π)^(-1)*B*(δ + A1.R - A2.R)))
					AtomNNs[n] = AtomNN(A1, A2, δ, N)
					n += 1	
				end
			end
		end
	end
	return AtomNNs
end

