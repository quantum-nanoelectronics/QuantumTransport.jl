using LinearAlgebra
include("Materials.jl")

# Takes in the parameters p, index vector (ix,iy,iz,isite (in unit cell), and iorb)
# returns the site-index in the full hamiltonian
function xyztoi(p,ivec, N::Vector{Int} = [0;0;0]) 
	# indexing 0 to N-1
	# in case the A₂ lattice vector != C*[0;1;0]
	diy = Int(round(p.SLa₂[1]/p.a₁[1]))*N[2]
	ix = mod(ivec[1],p.nx); iy = mod(ivec[2] + diy,p.ny); iz = mod(ivec[3],p.nz); isite = ivec[4]; iorb = ivec[5]
	return iorb + p.norb*isite + p.nsite*p.norb*ix + p.nsite*p.norb*p.nx*iy + p.nsite*p.norb*p.nx*p.ny*iz + 1
end

# Same as above, except returns the corresponding atomic position of each index vector 
# useful for calculating ∫A⋅δR peierls phase
function xyztor(p,ivec)
	ix = ivec[1]; iy = ivec[2]; iz = ivec[3]; isite = ivec[4];
	δdict = Dict(0 => p.A*[0.0; 0.0; 0.0], #In 1 
		     1 => p.A*[0.5; 0.5; 0.5]) #In 2
		     #2 => p.A*[0.0; 0.5; 0.8975-0.5], #Bi 1
		     #3 => p.A*[0.5; 0.0; 1.10248-0.5]) #Bi 2
	δ = δdict[isite]    
	R = p.a₁*ix + p.a₂*iy + p.a₃*iz + δ
	return R
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
			# at site 1->2 would be a_2,1
			# with site ⊗ spin, would be 
			H[2*NN.b-1, 2*NN.a-1] += NN.t[1,1]
			H[2*NN.b  , 2*NN.a-1] += NN.t[2,1]
			H[2*NN.b-1, 2*NN.a  ] += NN.t[1,2]
			H[2*NN.b  , 2*NN.a  ] += NN.t[2,2]
		end
	end
	return H, edgeNNs
end

rot(θ) = [cos(θ) -sin(θ); sin(θ) cos(θ)]


function genNNs(p) # all of the terms in the hamiltonian get added here, get back the relevant bonds
	NNs = Hopping[]
	function nextsite(isite::Int)
		return -2*isite + 1
	end
	# loop over each unit cell site in the superlattice
	for iy = 0:(p.ny-1)
		for iz = 0:(p.nz-1)
			for ix = 0:(p.nx-1)
				isite = 0
				for iorb = 0:(p.norb-1)
					ia = copy([ix,iy,iz,isite,iorb]);
					hoppingType! = hoppingDict[p.deviceMaterial]
					hoppingType!(p,NNs,ia)
				end
			end
		end
	end
	# now fix the designation for the vectors that hop out of the lattice
	# Will connect them around later using bloch's theorem to generate H(k) function
	for NN in NNs
		ib = [NN.ib[1],NN.ib[2],NN.ib[3]]
		# Δ(ib,ib reflected back into 1st lattice)
		pib = ib - [mod(ib[1],p.nx),mod(ib[2],p.ny),mod(ib[3],p.nz)]
		if(pib⋅pib != 0) # if vector is distinctly outside of 1st lattice
			NN.N = Int.([round(pib[1]/(p.nx)),round(pib[2]/p.ny),round(pib[3]/p.nz)])
			NN.b = xyztoi(p,NN.ib, NN.N)
			NN.edge = true
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
	return NNs
end